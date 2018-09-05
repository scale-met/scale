!-------------------------------------------------------------------------------
!> module COMMUNICATION
!!
!! @par Description
!!          MPI Communication module for Icosahedral A-grid
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_comm_icoA
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: COMM_setup
  public :: COMM_data_transfer
  public :: COMM_data_transfer_nopl
  public :: COMM_var

  interface COMM_data_transfer
     module procedure COMM_data_transfer_SP
     module procedure COMM_data_transfer_DP
  end interface COMM_data_transfer

  interface COMM_var
     module procedure COMM_var_SP
     module procedure COMM_var_DP
  end interface COMM_var

  public :: COMM_Stat_sum
  public :: COMM_Stat_sum_eachlayer
  public :: COMM_Stat_avg
  public :: COMM_Stat_max
  public :: COMM_Stat_min

  interface COMM_Stat_sum
     module procedure COMM_Stat_sum_SP
     module procedure COMM_Stat_sum_DP
  end interface COMM_Stat_sum

  interface COMM_Stat_sum_eachlayer
     module procedure COMM_Stat_sum_eachlayer_SP
     module procedure COMM_Stat_sum_eachlayer_DP
  end interface COMM_Stat_sum_eachlayer

  interface COMM_Stat_avg
     module procedure COMM_Stat_avg_SP
     module procedure COMM_Stat_avg_DP
  end interface COMM_Stat_avg

  interface COMM_Stat_max
     module procedure COMM_Stat_max_SP
     module procedure COMM_Stat_max_DP
  end interface COMM_Stat_max

  interface COMM_Stat_min
     module procedure COMM_Stat_min_SP
     module procedure COMM_Stat_min_DP
  end interface COMM_Stat_min

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,  public              :: COMM_pl = .true.

  integer,  public              :: COMM_datatype

  ! node-node communication/copy info
  integer,  public, allocatable :: Copy_info_r2r(:)   !> node-to-node information (intra-node copy)
  integer,  public, allocatable :: Recv_info_r2r(:,:) !> node-to-node information (receive)
  integer,  public, allocatable :: Send_info_r2r(:,:) !> node-to-node information (send)

  integer,  public, allocatable :: Copy_info_p2r(:)   !> node-to-node information (intra-node copy)
  integer,  public, allocatable :: Recv_info_p2r(:,:) !> node-to-node information (receive)
  integer,  public, allocatable :: Send_info_p2r(:,:) !> node-to-node information (send)

  integer,  public, allocatable :: Copy_info_r2p(:)   !> node-to-node information (intra-node copy)
  integer,  public, allocatable :: Recv_info_r2p(:,:) !> node-to-node information (receive)
  integer,  public, allocatable :: Send_info_r2p(:,:) !> node-to-node information (send)

  integer,  public, allocatable :: Singular_info(:)   !> node-to-node information (singular point)

  ! node-node communication/copy list
  integer,  public, allocatable :: Copy_list_r2r(:,:)   !> grid,local-region list (copy)
  integer,  public, allocatable :: Recv_list_r2r(:,:,:) !> grid,local-region list (receive)
  integer,  public, allocatable :: Send_list_r2r(:,:,:) !> grid,local-region list (send)

  integer,  public, allocatable :: Copy_list_p2r(:,:)   !> grid,local-region list (copy)
  integer,  public, allocatable :: Recv_list_p2r(:,:,:) !> grid,local-region list (receive)
  integer,  public, allocatable :: Send_list_p2r(:,:,:) !> grid,local-region list (send)

  integer,  public, allocatable :: Copy_list_r2p(:,:)   !> grid,local-region list (copy)
  integer,  public, allocatable :: Recv_list_r2p(:,:,:) !> grid,local-region list (receive)
  integer,  public, allocatable :: Send_list_r2p(:,:,:) !> grid,local-region list (send)

  integer,  public, allocatable :: Singular_list(:,:)   !> grid,local-region list (singular point)

  ! working buffer
  integer,  public, allocatable :: REQ_list(:)

  real(SP), public, allocatable :: sendbuf_r2r_SP(:,:)
  real(SP), public, allocatable :: recvbuf_r2r_SP(:,:)
  real(SP), public, allocatable :: sendbuf_p2r_SP(:,:)
  real(SP), public, allocatable :: recvbuf_p2r_SP(:,:)
  real(SP), public, allocatable :: sendbuf_r2p_SP(:,:)
  real(SP), public, allocatable :: recvbuf_r2p_SP(:,:)

  real(DP), public, allocatable :: sendbuf_r2r_DP(:,:)
  real(DP), public, allocatable :: recvbuf_r2r_DP(:,:)
  real(DP), public, allocatable :: sendbuf_p2r_DP(:,:)
  real(DP), public, allocatable :: recvbuf_p2r_DP(:,:)
  real(DP), public, allocatable :: sendbuf_r2p_DP(:,:)
  real(DP), public, allocatable :: recvbuf_r2p_DP(:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: COMM_list_generate
  private :: COMM_sortdest
  private :: COMM_sortdest_pl
  private :: COMM_sortdest_singular

  private :: COMM_debugtest

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private :: COMM_apply_barrier = .false. !< barrier option
  integer,  private :: COMM_varmax        = 50      !< maximum number of 3D variable for the communication buffer
  logical,  private :: debug              = .false. !< debug option
  logical,  private :: testonly           = .false. !< test  option

  ! send/recv relationship
  integer,  private, parameter   :: rellist_vindex = 6
  integer,  private, parameter   :: I_recv_grid = 1
  integer,  private, parameter   :: I_recv_rgn  = 2
  integer,  private, parameter   :: I_recv_prc  = 3
  integer,  private, parameter   :: I_send_grid = 4
  integer,  private, parameter   :: I_send_rgn  = 5
  integer,  private, parameter   :: I_send_prc  = 6

  integer,  private, allocatable :: rellist(:,:) !> send/recv relationship list
  integer,  private              :: rellist_nmax !> number of list

  ! node-node communication/copy info
  integer,  private, parameter   :: Recv_nlim = 20 !> number limit of rank to receive data
  integer,  private, parameter   :: Send_nlim = 20 !> number limit of rank to send    data

  integer,  private              :: Copy_nmax_r2r = 0
  integer,  private              :: Recv_nmax_r2r = 0
  integer,  private              :: Send_nmax_r2r = 0

  integer,  private              :: Copy_nmax_p2r = 0
  integer,  private              :: Recv_nmax_p2r = 0
  integer,  private              :: Send_nmax_p2r = 0

  integer,  private              :: Copy_nmax_r2p = 0
  integer,  private              :: Recv_nmax_r2p = 0
  integer,  private              :: Send_nmax_r2p = 0

  integer,  private              :: Singular_nmax = 0

  integer,  private, parameter   :: info_vindex = 3
  integer,  private, parameter   :: I_size     = 1
  integer,  private, parameter   :: I_prc_from = 2
  integer,  private, parameter   :: I_prc_to   = 3

  ! node-node communication/copy list
  integer,  private, parameter   :: list_vindex = 4
  integer,  private, parameter   :: I_grid_from = 1
  integer,  private, parameter   :: I_l_from    = 2
  integer,  private, parameter   :: I_grid_to   = 3
  integer,  private, parameter   :: I_l_to      = 4

  ! working buffer
  integer,  private              :: REQ_count

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine COMM_setup
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_RGN_r2p_pl, &
       I_NPL,          &
       I_SPL
    implicit none

    namelist / PARAM_COMM_ICOA / &
       COMM_apply_barrier, &
       COMM_varmax,        &
       debug,              &
       testonly

    integer  :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[comm]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COMM_ICOA,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** PARAM_COMM_ICOA is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist PARAM_COMM_ICOA. STOP.'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=PARAM_COMM_ICOA)

    if ( RP == DP ) then
       COMM_datatype = MPI_DOUBLE_PRECISION
    elseif( RP == SP ) then
       COMM_datatype = MPI_REAL
    else
       write(*,*) 'xxx precision is not supportd'
       call PRC_abort
    endif

    if (       PRC_RGN_r2p_pl(I_NPL) < 0 &
         .AND. PRC_RGN_r2p_pl(I_SPL) < 0 ) then
       COMM_pl = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '====== communication information ======'

    call COMM_list_generate

    call COMM_sortdest
    call COMM_sortdest_pl
    call COMM_sortdest_singular

    allocate( REQ_list( Recv_nmax_r2r + Send_nmax_r2r &
                      + Recv_nmax_p2r + Send_nmax_p2r &
                      + Recv_nmax_r2p + Send_nmax_r2p ) )

    if( testonly ) call COMM_debugtest

    return
  end subroutine COMM_setup

  !-----------------------------------------------------------------------------
  !> Generate inner grid -> halo communication list
  subroutine COMM_list_generate
    use scale_prc, only: &
       PRC_myrank
    use scale_prc_icoA, only: &
       PRC_RGN_local,    &
       PRC_RGN_r2lp,     &
       PRC_RGN_l2r,      &
       PRC_RGN_edge_tab, &
       PRC_RGN_vert_num, &
       PRC_RGN_vert_tab, &
       I_prc,            &
       I_RGNID,          &
       I_DIR,            &
       I_SW,             &
       I_NW,             &
       I_NE,             &
       I_SE,             &
       I_W,              &
       I_N,              &
       I_E,              &
       I_S
    implicit none

    integer  :: ginnar

    integer  :: prc, prc_rmt
    integer  :: rgnid, rgnid_rmt
    integer  :: i, j, i_rmt, j_rmt

    integer  :: n, l, cnt
    !---------------------------------------------------------------------------

    ginnar = ADM_gmax - ADM_gmin + 1

    allocate( rellist(rellist_vindex,ADM_gall*ADM_lall) )

    cnt = 0
    do l = 1, PRC_RGN_local
       rgnid = PRC_RGN_l2r(l)
       prc   = PRC_myrank

       !---< South West >---

       ! NE -> SW halo
       if ( PRC_RGN_edge_tab(I_DIR,I_SW,rgnid) == I_NE ) then
          rgnid_rmt = PRC_RGN_edge_tab(I_RGNID,I_SW,rgnid)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          do n = 1, ginnar
             cnt = cnt + 1

             i     = ADM_gmin - 1 + n
             j     = ADM_gmin - 1
             i_rmt = ADM_gmin - 1 + n
             j_rmt = ADM_gmax

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          enddo
       endif

       ! SE -> SW halo (Southern Hemisphere, Edge of diamond)
       if ( PRC_RGN_edge_tab(I_DIR,I_SW,rgnid) == I_SE ) then
          rgnid_rmt = PRC_RGN_edge_tab(I_RGNID,I_SW,rgnid)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          do n = 1, ginnar
             cnt = cnt + 1

             i     = ADM_gmin - 1 + n
             j     = ADM_gmin - 1
             i_rmt = ADM_gmax
             j_rmt = ADM_gmax + 1 - n ! reverse order

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          enddo
       endif

       !---< North West >---

       ! SE -> NW
       if ( PRC_RGN_edge_tab(I_DIR,I_NW,rgnid) == I_SE ) then
          rgnid_rmt = PRC_RGN_edge_tab(I_RGNID,I_NW,rgnid)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          do n = 1, ginnar
             cnt = cnt + 1

             i     = ADM_gmin - 1
             j     = ADM_gmin - 1 + n
             i_rmt = ADM_gmax
             j_rmt = ADM_gmin - 1 + n

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          enddo
       endif

       ! NE -> NW  (Northern Hemisphere, Edge of diamond)
       if ( PRC_RGN_edge_tab(I_DIR,I_NW,rgnid) == I_NE ) then
          rgnid_rmt = PRC_RGN_edge_tab(I_RGNID,I_NW,rgnid)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          do n = 1, ginnar
             cnt = cnt + 1

             i     = ADM_gmin - 1
             j     = ADM_gmin - 1 + n
             i_rmt = ADM_gmax + 1 - n ! reverse order
             j_rmt = ADM_gmax

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          enddo
       endif

       !---< North East >---

       ! SW -> NE
       if ( PRC_RGN_edge_tab(I_DIR,I_NE,rgnid) == I_SW ) then
          rgnid_rmt = PRC_RGN_edge_tab(I_RGNID,I_NE,rgnid)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          do n = 1, ginnar
             cnt = cnt + 1

             i     = ADM_gmin - 1 + n
             j     = ADM_gmax + 1
             i_rmt = ADM_gmin - 1 + n
             j_rmt = ADM_gmin

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          enddo
       endif

       ! NW -> NE  (Northern Hemisphere, Edge of diamond)
       if ( PRC_RGN_edge_tab(I_DIR,I_NE,rgnid) == I_NW ) then
          rgnid_rmt = PRC_RGN_edge_tab(I_RGNID,I_NE,rgnid)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          do n = 1, ginnar
             cnt = cnt + 1

             i     = ADM_gmin     + n ! shift 1 grid
             j     = ADM_gmax + 1
             i_rmt = ADM_gmin
             j_rmt = ADM_gmax + 1 - n ! reverse order

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          enddo
       endif

       !---< South East >---

       ! NW -> SE
       if ( PRC_RGN_edge_tab(I_DIR,I_SE,rgnid) == I_NW ) then
          rgnid_rmt = PRC_RGN_edge_tab(I_RGNID,I_SE,rgnid)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          do n = 1, ginnar
             cnt = cnt + 1

             i     = ADM_gmax + 1
             j     = ADM_gmin - 1 + n
             i_rmt = ADM_gmin
             j_rmt = ADM_gmin - 1 + n

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          enddo
       endif

       ! SW -> SE  (Southern Hemisphere, Edge of diamond)
       if ( PRC_RGN_edge_tab(I_DIR,I_SE,rgnid) == I_SW ) then
          rgnid_rmt = PRC_RGN_edge_tab(I_RGNID,I_SE,rgnid)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          do n = 1, ginnar
             cnt = cnt + 1

             i     = ADM_gmax + 1
             j     = ADM_gmin     + n ! shift 1 grid
             i_rmt = ADM_gmax + 1 - n ! reverse order
             j_rmt = ADM_gmin

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          enddo
       endif

       !---< Vertex : link to the next next region >---

       ! West Vertex
       if ( PRC_RGN_vert_num(I_W,rgnid) == 4 ) then
          rgnid_rmt = PRC_RGN_vert_tab(I_RGNID,I_W,rgnid,3)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          cnt = cnt + 1

          i = ADM_gmin - 1
          j = ADM_gmin - 1

          if ( PRC_RGN_vert_tab(I_DIR,I_W,rgnid,3) == I_N ) then
             i_rmt = ADM_gmin
             j_rmt = ADM_gmax
          elseif( PRC_RGN_vert_tab(I_DIR,I_W,rgnid,3) == I_E ) then
             i_rmt = ADM_gmax
             j_rmt = ADM_gmax
          elseif( PRC_RGN_vert_tab(I_DIR,I_W,rgnid,3) == I_S ) then
             i_rmt = ADM_gmax
             j_rmt = ADM_gmin
          endif

          rellist(I_recv_grid,cnt) = suf(i,j)
          rellist(I_recv_rgn, cnt) = rgnid
          rellist(I_recv_prc, cnt) = prc
          rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
          rellist(I_send_rgn, cnt) = rgnid_rmt
          rellist(I_send_prc, cnt) = prc_rmt
       endif

       ! North Vertex
       if ( PRC_RGN_vert_num(I_N,rgnid) == 4 ) then
          rgnid_rmt = PRC_RGN_vert_tab(I_RGNID,I_N,rgnid,3)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          ! known as north pole point
          if ( PRC_RGN_vert_tab(I_DIR,I_N,rgnid,3) == I_W ) then
             cnt = cnt + 1

             i     = ADM_gmin
             j     = ADM_gmax + 1
             i_rmt = ADM_gmin
             j_rmt = ADM_gmin

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          endif

          ! unused vertex point
          cnt = cnt + 1

          i = ADM_gmin - 1
          j = ADM_gmax + 1

          if ( PRC_RGN_vert_tab(I_DIR,I_N,rgnid,3) == I_W ) then
             i_rmt = ADM_gmin
             j_rmt = ADM_gmin + 1
          elseif( PRC_RGN_vert_tab(I_DIR,I_N,rgnid,3) == I_N ) then
             i_rmt = ADM_gmin
             j_rmt = ADM_gmax
          elseif( PRC_RGN_vert_tab(I_DIR,I_N,rgnid,3) == I_E ) then
             i_rmt = ADM_gmax
             j_rmt = ADM_gmax
          elseif( PRC_RGN_vert_tab(I_DIR,I_N,rgnid,3) == I_S ) then
             i_rmt = ADM_gmax
             j_rmt = ADM_gmin
          endif

          rellist(I_recv_grid,cnt) = suf(i,j)
          rellist(I_recv_rgn, cnt) = rgnid
          rellist(I_recv_prc, cnt) = prc
          rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
          rellist(I_send_rgn, cnt) = rgnid_rmt
          rellist(I_send_prc, cnt) = prc_rmt
       endif

       ! East Vertex
       if ( PRC_RGN_vert_num(I_E,rgnid) == 4 ) then
          if ( PRC_RGN_vert_tab(I_DIR,I_E,rgnid,3) == I_W ) then
             rgnid_rmt = PRC_RGN_vert_tab(I_RGNID,I_E,rgnid,3)
             prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

             cnt = cnt + 1

             i     = ADM_gmax + 1
             j     = ADM_gmax + 1
             i_rmt = ADM_gmin
             j_rmt = ADM_gmin

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          endif
       endif

       ! South Vertex
       if ( PRC_RGN_vert_num(I_S,rgnid) == 4 ) then
          rgnid_rmt = PRC_RGN_vert_tab(I_RGNID,I_S,rgnid,3)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          ! known as south pole point
          if ( PRC_RGN_vert_tab(I_DIR,I_S,rgnid,3) == I_W ) then
             cnt = cnt + 1

             i     = ADM_gmax + 1
             j     = ADM_gmin
             i_rmt = ADM_gmin
             j_rmt = ADM_gmin

             rellist(I_recv_grid,cnt) = suf(i,j)
             rellist(I_recv_rgn, cnt) = rgnid
             rellist(I_recv_prc, cnt) = prc
             rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
             rellist(I_send_rgn, cnt) = rgnid_rmt
             rellist(I_send_prc, cnt) = prc_rmt
          endif

          ! unused vertex point
          cnt = cnt + 1

          i = ADM_gmax + 1
          j = ADM_gmin - 1

          if ( PRC_RGN_vert_tab(I_DIR,I_S,rgnid,3) == I_W ) then
             i_rmt = ADM_gmin + 1
             j_rmt = ADM_gmin
          elseif( PRC_RGN_vert_tab(I_DIR,I_S,rgnid,3) == I_N ) then
             i_rmt = ADM_gmin
             j_rmt = ADM_gmax
          elseif( PRC_RGN_vert_tab(I_DIR,I_S,rgnid,3) == I_E ) then
             i_rmt = ADM_gmax
             j_rmt = ADM_gmax
          elseif( PRC_RGN_vert_tab(I_DIR,I_S,rgnid,3) == I_S ) then
             i_rmt = ADM_gmax
             j_rmt = ADM_gmin
          endif

          rellist(I_recv_grid,cnt) = suf(i,j)
          rellist(I_recv_rgn, cnt) = rgnid
          rellist(I_recv_prc, cnt) = prc
          rellist(I_send_grid,cnt) = suf(i_rmt,j_rmt)
          rellist(I_send_rgn, cnt) = rgnid_rmt
          rellist(I_send_prc, cnt) = prc_rmt
       endif
    enddo ! loop l

    rellist_nmax = cnt

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** rellist_nmax:', rellist_nmax

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*) '--- Relation Table'
       if( IO_L ) write(IO_FID_LOG,'(7(A10))') 'Count', '|recv_grid', '| recv_rgn', '| recv_prc', &
                                                        '|send_grid', '| send_rgn', '| send_prc'
       do cnt = 1, rellist_nmax
          if( IO_L ) write(IO_FID_LOG,'(7(I10))') cnt, rellist(:,cnt)
       enddo
    endif

    return
  end subroutine COMM_list_generate

  !-----------------------------------------------------------------------------
  !> Sort data destination for region <-> region
  subroutine COMM_sortdest
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_myrank,           &
       PRC_nprocs
    use scale_prc_icoA, only: &
       PRC_RGN_lp2r, &
       PRC_RGN_r2lp, &
       I_l
    implicit none

    integer  :: sendbuf1(1)
    integer  :: recvbuf1(1)

    integer, allocatable :: sendbuf_info(:)
    integer, allocatable :: recvbuf_info(:)
    integer, allocatable :: sendbuf_list(:,:,:)
    integer, allocatable :: recvbuf_list(:,:,:)
    integer, allocatable :: REQ_list_r2r(:)

    integer  :: Recv_nglobal_r2r
    integer  :: Send_size_nglobal

    integer  :: cnt, irank, ipos
    integer  :: totalsize, rank, tag

    integer  :: i_from, j_from, r_from, g_from, l_from, p_from
    integer  :: i_to, j_to, r_to, g_to, l_to, p_to

    integer  :: ierr
    integer  :: n, p
    !---------------------------------------------------------------------------

    allocate( Copy_info_r2r(info_vindex) )
    allocate( Recv_info_r2r(info_vindex,Recv_nlim) )
    allocate( Send_info_r2r(info_vindex,Send_nlim) )
    Copy_info_r2r(:)        = -1
    Recv_info_r2r(:,:)      = -1
    Send_info_r2r(:,:)      = -1
    Copy_info_r2r(I_size)   = 0
    Recv_info_r2r(I_size,:) = 0
    Send_info_r2r(I_size,:) = 0

    allocate( Copy_list_r2r(list_vindex,rellist_nmax) )
    allocate( Recv_list_r2r(list_vindex,rellist_nmax,Recv_nlim) )
    allocate( Send_list_r2r(list_vindex,rellist_nmax,Send_nlim) )
    Copy_list_r2r(:,:)   = -1
    Recv_list_r2r(:,:,:) = -1
    Send_list_r2r(:,:,:) = -1

    ! sorting list according to destination
    do cnt = 1, rellist_nmax

       if ( rellist(I_recv_prc,cnt) == rellist(I_send_prc,cnt) ) then ! no communication
          Copy_info_r2r(I_size) = Copy_info_r2r(I_size) + 1
          ipos                  = Copy_info_r2r(I_size)

          Copy_list_r2r(I_grid_from,ipos) = rellist(I_send_grid,cnt)
          Copy_list_r2r(I_l_from   ,ipos) = PRC_RGN_r2lp(I_l,rellist(I_send_rgn,cnt))
          Copy_list_r2r(I_grid_to  ,ipos) = rellist(I_recv_grid,cnt)
          Copy_list_r2r(I_l_to     ,ipos) = PRC_RGN_r2lp(I_l,rellist(I_recv_rgn,cnt))
       else ! node-to-node communication
          !--- search existing rank id (identify by prc_from)
          irank = -1
          do n = 1, Recv_nmax_r2r
             if ( Recv_info_r2r(I_prc_from,n) == rellist(I_send_prc,cnt) ) then
                irank = n
                exit
             endif
          enddo

          if ( irank < 0 ) then ! register new rank id
             Recv_nmax_r2r = Recv_nmax_r2r + 1
             irank         = Recv_nmax_r2r

             Recv_info_r2r(I_prc_from,irank) = rellist(I_send_prc,cnt)
             Recv_info_r2r(I_prc_to  ,irank) = rellist(I_recv_prc,cnt)
          endif

          Recv_info_r2r(I_size,irank) = Recv_info_r2r(I_size,irank) + 1
          ipos                        = Recv_info_r2r(I_size,irank)

          Recv_list_r2r(I_grid_from,ipos,irank) = rellist(I_send_grid,cnt)
          Recv_list_r2r(I_l_from   ,ipos,irank) = PRC_RGN_r2lp(I_l,rellist(I_send_rgn,cnt))
          Recv_list_r2r(I_grid_to  ,ipos,irank) = rellist(I_recv_grid,cnt)
          Recv_list_r2r(I_l_to     ,ipos,irank) = PRC_RGN_r2lp(I_l,rellist(I_recv_rgn,cnt))
       endif
    enddo

    if ( Copy_info_r2r(I_size) > 0 ) then
       Copy_nmax_r2r = 1
       Copy_info_r2r(I_prc_from) = PRC_myrank
       Copy_info_r2r(I_prc_to  ) = PRC_myrank
    endif



    ! get maximum number of rank to communication
    sendbuf1(1) = Recv_nmax_r2r

    call MPI_Allreduce( sendbuf1(1),          & ! [IN]
                        recvbuf1(1),          & ! [OUT]
                        1,                    & ! [IN]
                        MPI_INTEGER,          & ! [IN]
                        MPI_MAX,              & ! [IN]
                        PRC_LOCAL_COMM_WORLD, & ! [IN]
                        ierr                  ) ! [OUT]

    Recv_nglobal_r2r = recvbuf1(1)

    ! Distribute receive request from each rank to all
    allocate( sendbuf_info(info_vindex*Recv_nglobal_r2r) )
    allocate( recvbuf_info(info_vindex*Recv_nglobal_r2r*PRC_nprocs) )
    sendbuf_info(:) = -1

    do irank = 1, Recv_nmax_r2r
       n = (irank-1) * info_vindex

       sendbuf_info(n+I_size    ) = Recv_info_r2r(I_size    ,irank)
       sendbuf_info(n+I_prc_from) = Recv_info_r2r(I_prc_from,irank)
       sendbuf_info(n+I_prc_to  ) = Recv_info_r2r(I_prc_to  ,irank)
    enddo

    totalsize = info_vindex * Recv_nglobal_r2r

    if ( totalsize > 0 ) then
       call MPI_Allgather( sendbuf_info(1),      & ! [IN]
                           totalsize,            & ! [IN]
                           MPI_INTEGER,          & ! [IN]
                           recvbuf_info(1),      & ! [OUT]
                           totalsize,            & ! [IN]
                           MPI_INTEGER,          & ! [IN]
                           PRC_LOCAL_COMM_WORLD, & ! [IN]
                           ierr                  ) ! [OUT]
    endif

    Send_size_nglobal = 0

    ! Accept receive request to my rank
    do p  = 1, Recv_nglobal_r2r*PRC_nprocs
       n = (p-1) * info_vindex

       if ( recvbuf_info(n+I_prc_from) == PRC_myrank ) then
          Send_nmax_r2r = Send_nmax_r2r + 1
          irank         = Send_nmax_r2r

          Send_info_r2r(I_size    ,irank) = recvbuf_info(n+I_size    )
          Send_info_r2r(I_prc_from,irank) = recvbuf_info(n+I_prc_from)
          Send_info_r2r(I_prc_to  ,irank) = recvbuf_info(n+I_prc_to  )
       endif

       Send_size_nglobal = max( Send_size_nglobal, recvbuf_info(n+I_size) )
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Recv_nmax_r2r(global) = ', Recv_nglobal_r2r
    if( IO_L ) write(IO_FID_LOG,*) '*** Recv_nmax_r2r(local)  = ', Recv_nmax_r2r
    if( IO_L ) write(IO_FID_LOG,*) '*** Send_nmax_r2r(local)  = ', Send_nmax_r2r
    if( IO_L ) write(IO_FID_LOG,*) '*** Send_size_r2r(global) = ', Send_size_nglobal
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(A)')             '|---------------------------------------'
    if( IO_L ) write(IO_FID_LOG,'(A)')             '|               size  prc_from    prc_to'
    if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))')    '| Copy_r2r', Copy_info_r2r(:)
    do irank = 1, Recv_nmax_r2r
       if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))') '| Recv_r2r', Recv_info_r2r(:,irank)
    enddo
    do irank = 1, Send_nmax_r2r
       if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))') '| Send_r2r', Send_info_r2r(:,irank)
    enddo

    ! Communicate detailed information in each pair
    allocate( REQ_list_r2r(Recv_nmax_r2r+Send_nmax_r2r)  )

    allocate( sendbuf_list(list_vindex,Send_size_nglobal,Recv_nmax_r2r) )
    allocate( recvbuf_list(list_vindex,Send_size_nglobal,Send_nmax_r2r) )
    sendbuf_list(:,:,:) = -1

    REQ_count = 0

    do irank = 1, Send_nmax_r2r
       REQ_count = REQ_count + 1
       totalsize = Send_info_r2r(I_size,irank) * list_vindex
       rank      = Send_info_r2r(I_prc_to  ,irank)
       tag       = Send_info_r2r(I_prc_from,irank)

       call MPI_IRECV( recvbuf_list(1,1,irank), & ! [OUT]
                       totalsize,               & ! [IN]
                       MPI_INTEGER,             & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list_r2r(REQ_count), & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    do irank = 1, Recv_nmax_r2r
       do ipos = 1, Recv_info_r2r(I_size,irank)
          sendbuf_list(:,ipos,irank) = Recv_list_r2r(:,ipos,irank)
       enddo

       REQ_count = REQ_count + 1
       totalsize = Recv_info_r2r(I_size,irank) * list_vindex
       rank      = Recv_info_r2r(I_prc_from,irank)
       tag       = Recv_info_r2r(I_prc_from,irank)

       call MPI_ISEND( sendbuf_list(1,1,irank), & ! [IN]
                       totalsize,               & ! [IN]
                       MPI_INTEGER,             & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list_r2r(REQ_count), & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- wait packets
    if ( Recv_nmax_r2r+Send_nmax_r2r > 0 ) then
       call MPI_WAITALL( Recv_nmax_r2r+Send_nmax_r2r, & ! [IN]
                         REQ_list_r2r(:),             & ! [IN]
                         MPI_STATUSES_IGNORE,         & ! [OUT]
                         ierr                         ) ! [OUT]
    endif

    do irank = 1, Send_nmax_r2r
       do ipos = 1, Send_info_r2r(I_size,irank)
          Send_list_r2r(:,ipos,irank) = recvbuf_list(:,ipos,irank)
       enddo
    enddo

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Copy_list_r2r'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(13(A6))') 'number', &
                                               '|ifrom','|jfrom','|rfrom','|gfrom','|lfrom','|pfrom', &
                                               '|  ito','|  jto','|  rto','|  gto','|  lto','|  pto'
       do ipos = 1, Copy_info_r2r(I_size)
          g_from = Copy_list_r2r(I_grid_from,ipos)
          l_from = Copy_list_r2r(I_l_from,ipos)
          p_from = Copy_info_r2r(I_prc_from)
          i_from = mod(g_from-1,ADM_gall_1d) + 1
          j_from = (g_from-i_from) / ADM_gall_1d + 1
          r_from = PRC_RGN_lp2r(l_from,p_from)
          g_to   = Copy_list_r2r(I_grid_to,ipos)
          l_to   = Copy_list_r2r(I_l_to,ipos)
          p_to   = Copy_info_r2r(I_prc_to)
          i_to   = mod(g_to-1,ADM_gall_1d) + 1
          j_to   = (g_to-i_to) / ADM_gall_1d + 1
          r_to   = PRC_RGN_lp2r(l_to,p_to)

          if( IO_L ) write(IO_FID_LOG,'(13(I6))') ipos, i_from, j_from, r_from, g_from, l_from, p_from, &
                                                        i_to  , j_to  , r_to  , g_to  , l_to  , p_to
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Recv_list_r2r'
       do irank = 1, Recv_nmax_r2r
          if( IO_L ) write(IO_FID_LOG,'(13(A6))') 'number', &
                                                  '|ifrom','|jfrom','|rfrom','|gfrom','|lfrom','|pfrom', &
                                                  '|  ito','|  jto','|  rto','|  gto','|  lto','|  pto'
          do ipos = 1, Recv_info_r2r(I_size,irank)
             g_from = Recv_list_r2r(I_grid_from,ipos,irank)
             l_from = Recv_list_r2r(I_l_from,ipos,irank)
             p_from = Recv_info_r2r(I_prc_from,irank)
             i_from = mod(g_from-1,ADM_gall_1d) + 1
             j_from = (g_from-i_from) / ADM_gall_1d + 1
             r_from = PRC_RGN_lp2r(l_from,p_from)
             g_to   = Recv_list_r2r(I_grid_to,ipos,irank)
             l_to   = Recv_list_r2r(I_l_to,ipos,irank)
             p_to   = Recv_info_r2r(I_prc_to,irank)
             i_to   = mod(g_to-1,ADM_gall_1d) + 1
             j_to   = (g_to-i_to) / ADM_gall_1d + 1
             r_to   = PRC_RGN_lp2r(l_to,p_to)

             if( IO_L ) write(IO_FID_LOG,'(13(I6))') ipos, i_from, j_from, r_from, g_from, l_from, p_from, &
                                                           i_to  , j_to  , r_to  , g_to  , l_to  , p_to
          enddo
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Send_list_r2r'
       do irank = 1, Send_nmax_r2r
          if( IO_L ) write(IO_FID_LOG,'(13(A6))') 'number', &
                                                  '|ifrom','|jfrom','|rfrom','|gfrom','|lfrom','|pfrom', &
                                                  '|  ito','|  jto','|  rto','|  gto','|  lto','|  pto'
          do ipos = 1, Send_info_r2r(I_size,irank)
             g_from = Send_list_r2r(I_grid_from,ipos,irank)
             l_from = Send_list_r2r(I_l_from,ipos,irank)
             p_from = Send_info_r2r(I_prc_from,irank)
             i_from = mod(g_from-1,ADM_gall_1d) + 1
             j_from = (g_from-i_from) / ADM_gall_1d + 1
             r_from = PRC_RGN_lp2r(l_from,p_from)
             g_to   = Send_list_r2r(I_grid_to,ipos,irank)
             l_to   = Send_list_r2r(I_l_to,ipos,irank)
             p_to   = Send_info_r2r(I_prc_to,irank)
             i_to   = mod(g_to-1,ADM_gall_1d) + 1
             j_to   = (g_to-i_to) / ADM_gall_1d + 1
             r_to   = PRC_RGN_lp2r(l_to,p_to)

             if( IO_L ) write(IO_FID_LOG,'(13(I6))') ipos, i_from, j_from, r_from, g_from, l_from, p_from, &
                                                           i_to  , j_to  , r_to  , g_to  , l_to  , p_to
          enddo
       enddo
    endif

    allocate( sendbuf_r2r_SP(Send_size_nglobal*ADM_kall*COMM_varmax,Send_nmax_r2r) )
    allocate( recvbuf_r2r_SP(Send_size_nglobal*ADM_kall*COMM_varmax,Recv_nmax_r2r) )
    allocate( sendbuf_r2r_DP(Send_size_nglobal*ADM_kall*COMM_varmax,Send_nmax_r2r) )
    allocate( recvbuf_r2r_DP(Send_size_nglobal*ADM_kall*COMM_varmax,Recv_nmax_r2r) )

    return
  end subroutine COMM_sortdest

  !-----------------------------------------------------------------------------
  !> Setup data communication for region <-> pole
  subroutine COMM_sortdest_pl
    use scale_prc, only: &
       PRC_myrank
    use scale_prc_icoA, only: &
       PRC_have_pl,         &
       PRC_RGN_local,       &
       PRC_RGN_vert_num,    &
       PRC_RGN_vert_tab_pl, &
       PRC_RGN_lp2r,        &
       PRC_RGN_r2lp,        &
       PRC_RGN_l2r,         &
       PRC_RGN_r2p_pl,      &
       I_l,                 &
       I_prc,               &
       I_RGNID,             &
       I_N,                 &
       I_S,                 &
       I_NPL,               &
       I_SPL
    implicit none

    integer  :: prc, prc_rmt
    integer  :: rgnid, rgnid_rmt
    integer  :: check_vert_num
    integer  :: pl_to

    integer  :: irank, ipos

    integer, parameter :: Send_size_nglobal_pl = 10

    integer  :: l, l_pl, n, v, vv
    integer  :: i_from, j_from, r_from, g_from, l_from, p_from
    integer  :: i_to, j_to, r_to, g_to, l_to, p_to
    !---------------------------------------------------------------------------

    allocate( Copy_info_p2r(info_vindex) )
    allocate( Recv_info_p2r(info_vindex,Recv_nlim) )
    allocate( Send_info_p2r(info_vindex,Send_nlim) )
    Copy_info_p2r(:)        = -1
    Recv_info_p2r(:,:)      = -1
    Send_info_p2r(:,:)      = -1
    Copy_info_p2r(I_size)   = 0
    Recv_info_p2r(I_size,:) = 0
    Send_info_p2r(I_size,:) = 0

    allocate( Copy_list_p2r(list_vindex,Send_size_nglobal_pl) )
    allocate( Recv_list_p2r(list_vindex,Send_size_nglobal_pl,Recv_nlim) )
    allocate( Send_list_p2r(list_vindex,Send_size_nglobal_pl,Send_nlim) )
    Copy_list_p2r(:,:)   = -1
    Recv_list_p2r(:,:,:) = -1
    Send_list_p2r(:,:,:) = -1

    allocate( Copy_info_r2p(info_vindex) )
    allocate( Recv_info_r2p(info_vindex,Recv_nlim) )
    allocate( Send_info_r2p(info_vindex,Send_nlim) )
    Copy_info_r2p(:)        = -1
    Recv_info_r2p(:,:)      = -1
    Send_info_r2p(:,:)      = -1
    Copy_info_r2p(I_size)   = 0
    Recv_info_r2p(I_size,:) = 0
    Send_info_r2p(I_size,:) = 0

    allocate( Copy_list_r2p(list_vindex,Send_size_nglobal_pl) )
    allocate( Recv_list_r2p(list_vindex,Send_size_nglobal_pl,Recv_nlim) )
    allocate( Send_list_r2p(list_vindex,Send_size_nglobal_pl,Send_nlim) )
    Copy_list_r2p(:,:)   = -1
    Recv_list_r2p(:,:,:) = -1
    Send_list_r2p(:,:,:) = -1

    ! Search in regular region
    do l = 1, PRC_RGN_local
       rgnid = PRC_RGN_l2r(l)
       prc   = PRC_myrank

       do l_pl = I_NPL, I_SPL
          rgnid_rmt = l_pl
          prc_rmt   = PRC_RGN_r2p_pl(rgnid_rmt)

          if ( rgnid_rmt == I_NPL ) then
             check_vert_num = PRC_RGN_vert_num(I_N,rgnid)
             i_from = ADM_gmin
             j_from = ADM_gmax
             i_to   = ADM_gmin
             j_to   = ADM_gmax + 1
          elseif( rgnid_rmt == I_SPL ) then
             check_vert_num = PRC_RGN_vert_num(I_S,rgnid)
             i_from = ADM_gmax
             j_from = ADM_gmin
             i_to   = ADM_gmax + 1
             j_to   = ADM_gmin
          endif

          if ( check_vert_num == ADM_vlink ) then

             ! search destination in the pole halo
             do v = 1, ADM_vlink
                if ( rgnid == PRC_RGN_vert_tab_pl(I_RGNID,rgnid_rmt,v) ) then
                   pl_to = v
                   if( pl_to < ADM_gmin_pl ) pl_to = ADM_gmax_pl
                   exit
               endif
             enddo

             if ( prc == prc_rmt ) then ! no communication
                ! copy region inner -> pole halo
                Copy_info_r2p(I_size) = Copy_info_r2p(I_size) + 1
                ipos                  = Copy_info_r2p(I_size)

                Copy_list_r2p(I_grid_from,ipos) = suf(i_from,j_from)
                Copy_list_r2p(I_l_from   ,ipos) = l
                Copy_list_r2p(I_grid_to  ,ipos) = pl_to
                Copy_list_r2p(I_l_to     ,ipos) = l_pl
             else ! node-to-node communication
                ! receive pole center -> region halo
                irank = -1
                do n = 1, Recv_nmax_p2r
                   if ( Recv_info_p2r(I_prc_from,n) == prc_rmt ) then
                      irank = n
                      exit
                   endif
                enddo

                if ( irank < 0 ) then ! register new rank id
                   Recv_nmax_p2r = Recv_nmax_p2r + 1
                   irank         = Recv_nmax_p2r

                   Recv_info_p2r(I_prc_from,irank) = prc_rmt
                   Recv_info_p2r(I_prc_to  ,irank) = prc
                endif

                Recv_info_p2r(I_size,irank) = Recv_info_p2r(I_size,irank) + 1
                ipos                        = Recv_info_p2r(I_size,irank)

                Recv_list_p2r(I_grid_from,ipos,irank) = ADM_gslf_pl
                Recv_list_p2r(I_l_from   ,ipos,irank) = l_pl
                Recv_list_p2r(I_grid_to  ,ipos,irank) = suf(i_to,j_to)
                Recv_list_p2r(I_l_to     ,ipos,irank) = l

                ! send region inner -> pole halo
                irank = -1
                do n = 1, Send_nmax_r2p
                   if ( Send_info_r2p(I_prc_to,n) == prc_rmt ) then
                      irank = n
                      exit
                   endif
                enddo

                if ( irank < 0 ) then ! register new rank id
                   Send_nmax_r2p = Send_nmax_r2p + 1
                   irank         = Send_nmax_r2p

                   Send_info_r2p(I_prc_from,irank) = prc
                   Send_info_r2p(I_prc_to  ,irank) = prc_rmt
                endif

                Send_info_r2p(I_size,irank) = Send_info_r2p(I_size,irank) + 1
                ipos                        = Send_info_r2p(I_size,irank)

                Send_list_r2p(I_grid_from,ipos,irank) = suf(i_from,j_from)
                Send_list_r2p(I_l_from   ,ipos,irank) = l
                Send_list_r2p(I_grid_to  ,ipos,irank) = pl_to
                Send_list_r2p(I_l_to     ,ipos,irank) = l_pl
             endif
          endif

       enddo ! loop l_pl
    enddo ! loop l

    ! Search in pole
    if ( PRC_have_pl ) then

    do l_pl = I_NPL, I_SPL
       rgnid = l_pl
       prc   = PRC_myrank

       do v = ADM_vlink+1, 2, -1
          vv = v
          if( v == ADM_vlink+1 ) vv = 1

          rgnid_rmt = PRC_RGN_vert_tab_pl(I_RGNID,rgnid,vv)
          prc_rmt   = PRC_RGN_r2lp(I_prc,rgnid_rmt)

          if ( rgnid == I_NPL ) then
             i_from = ADM_gmin
             j_from = ADM_gmax
             i_to   = ADM_gmin
             j_to   = ADM_gmax + 1
          elseif( rgnid == I_SPL ) then
             i_from = ADM_gmax
             j_from = ADM_gmin
             i_to   = ADM_gmax + 1
             j_to   = ADM_gmin
          endif

          pl_to = vv
          if( pl_to < ADM_gmin_pl ) pl_to = ADM_gmax_pl

          if ( prc == prc_rmt ) then ! no communication
             ! copy pole center -> region halo
             Copy_info_p2r(I_size) = Copy_info_p2r(I_size) + 1
             ipos                  = Copy_info_p2r(I_size)

             Copy_list_p2r(I_grid_from,ipos) = ADM_gslf_pl
             Copy_list_p2r(I_l_from   ,ipos) = l_pl
             Copy_list_p2r(I_grid_to  ,ipos) = suf(i_to,j_to)
             Copy_list_p2r(I_l_to     ,ipos) = PRC_RGN_r2lp(I_l,rgnid_rmt)
          else ! node-to-node communication
             ! send region inner -> pole halo
             irank = -1
             do n = 1, Recv_nmax_r2p
                if ( Recv_info_r2p(I_prc_from,n) == prc_rmt ) then
                   irank = n
                   exit
                endif
             enddo

             if ( irank < 0 ) then ! register new rank id
                Recv_nmax_r2p = Recv_nmax_r2p + 1
                irank         = Recv_nmax_r2p

                Recv_info_r2p(I_prc_from,irank) = prc_rmt
                Recv_info_r2p(I_prc_to  ,irank) = prc
             endif

             Recv_info_r2p(I_size,irank) = Recv_info_r2p(I_size,irank) + 1
             ipos                        = Recv_info_r2p(I_size,irank)

             Recv_list_r2p(I_grid_from,ipos,irank) = suf(i_from,j_from)
             Recv_list_r2p(I_l_from   ,ipos,irank) = PRC_RGN_r2lp(I_l,rgnid_rmt)
             Recv_list_r2p(I_grid_to  ,ipos,irank) = pl_to
             Recv_list_r2p(I_l_to     ,ipos,irank) = l_pl

             ! receive pole center -> region halo
             irank = -1
             do n = 1, Send_nmax_p2r
                if ( Send_info_p2r(I_prc_to,n) == prc_rmt ) then
                   irank = n
                   exit
                endif
             enddo

             if ( irank < 0 ) then ! register new rank id
                Send_nmax_p2r = Send_nmax_p2r + 1
                irank         = Send_nmax_p2r

                Send_info_p2r(I_prc_from,irank) = prc
                Send_info_p2r(I_prc_to  ,irank) = prc_rmt
             endif

             Send_info_p2r(I_size,irank) = Send_info_p2r(I_size,irank) + 1
             ipos                        = Send_info_p2r(I_size,irank)

             Send_list_p2r(I_grid_from,ipos,irank) = ADM_gslf_pl
             Send_list_p2r(I_l_from   ,ipos,irank) = l_pl
             Send_list_p2r(I_grid_to  ,ipos,irank) = suf(i_to,j_to)
             Send_list_p2r(I_l_to     ,ipos,irank) = PRC_RGN_r2lp(I_l,rgnid_rmt)
          endif

       enddo ! loop pole halo

    enddo ! loop l_pl

    if ( Copy_info_p2r(I_size) > 0 ) then
       Copy_nmax_p2r = 1
       Copy_info_p2r(I_prc_from) = PRC_myrank
       Copy_info_p2r(I_prc_to  ) = PRC_myrank
    endif

    if ( Copy_info_r2p(I_size) > 0 ) then
       Copy_nmax_r2p = 1
       Copy_info_r2p(I_prc_from) = PRC_myrank
       Copy_info_r2p(I_prc_to  ) = PRC_myrank
    endif

    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Recv_nmax_p2r(local)  = ', Recv_nmax_p2r
    if( IO_L ) write(IO_FID_LOG,*) '*** Send_nmax_p2r(local)  = ', Send_nmax_p2r
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(A)')             '|---------------------------------------'
    if( IO_L ) write(IO_FID_LOG,'(A)')             '|               size  prc_from    prc_to'
    if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))')    '| Copy_p2r', Copy_info_p2r(:)
    do irank = 1, Recv_nmax_p2r
       if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))') '| Recv_p2r', Recv_info_p2r(:,irank)
    enddo
    do irank = 1, Send_nmax_p2r
       if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))') '| Send_p2r', Send_info_p2r(:,irank)
    enddo

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Copy_list_p2r'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(11(A6))') 'number', &
                                                                 '|rfrom','|gfrom','|lfrom','|pfrom', &
                                               '|  ito','|  jto','|  rto','|  gto','|  lto','|  pto'
       do ipos = 1, Copy_info_p2r(I_size)
          g_from = Copy_list_p2r(I_grid_from,ipos)
          l_from = Copy_list_p2r(I_l_from,ipos)
          p_from = Copy_info_p2r(I_prc_from)
          r_from = l_from
          g_to   = Copy_list_p2r(I_grid_to,ipos)
          l_to   = Copy_list_p2r(I_l_to,ipos)
          p_to   = Copy_info_p2r(I_prc_to)
          i_to   = mod(g_to-1,ADM_gall_1d) + 1
          j_to   = (g_to-i_to) / ADM_gall_1d + 1
          r_to   = PRC_RGN_lp2r(l_to,p_to)

          if( IO_L ) write(IO_FID_LOG,'(11(I6))') ipos,                 r_from, g_from, l_from, p_from, &
                                                        i_to  , j_to  , r_to  , g_to  , l_to  , p_to
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Recv_list_p2r'
       do irank = 1, Recv_nmax_p2r
          if( IO_L ) write(IO_FID_LOG,'(11(A6))') 'number', &
                                                                    '|rfrom','|gfrom','|lfrom','|pfrom', &
                                                  '|  ito','|  jto','|  rto','|  gto','|  lto','|  pto'
          do ipos = 1, Recv_info_p2r(I_size,irank)
             g_from = Recv_list_p2r(I_grid_from,ipos,irank)
             l_from = Recv_list_p2r(I_l_from,ipos,irank)
             p_from = Recv_info_p2r(I_prc_from,irank)
             r_from = l_from
             g_to   = Recv_list_p2r(I_grid_to,ipos,irank)
             l_to   = Recv_list_p2r(I_l_to,ipos,irank)
             p_to   = Recv_info_p2r(I_prc_to,irank)
             i_to   = mod(g_to-1,ADM_gall_1d) + 1
             j_to   = (g_to-i_to) / ADM_gall_1d + 1
             r_to   = PRC_RGN_lp2r(l_to,p_to)

             if( IO_L ) write(IO_FID_LOG,'(11(I6))') ipos,                 r_from, g_from, l_from, p_from, &
                                                           i_to  , j_to  , r_to  , g_to  , l_to  , p_to
          enddo
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Send_list_p2r'
       do irank = 1, Send_nmax_p2r
          if( IO_L ) write(IO_FID_LOG,'(11(A6))') 'number', &
                                                                    '|rfrom','|gfrom','|lfrom','|pfrom', &
                                                  '|  ito','|  jto','|  rto','|  gto','|  lto','|  pto'
          do ipos = 1, Send_info_p2r(I_size,irank)
             g_from = Send_list_p2r(I_grid_from,ipos,irank)
             l_from = Send_list_p2r(I_l_from,ipos,irank)
             p_from = Send_info_p2r(I_prc_from,irank)
             r_from = l_from
             g_to   = Send_list_p2r(I_grid_to,ipos,irank)
             l_to   = Send_list_p2r(I_l_to,ipos,irank)
             p_to   = Send_info_p2r(I_prc_to,irank)
             i_to   = mod(g_to-1,ADM_gall_1d) + 1
             j_to   = (g_to-i_to) / ADM_gall_1d + 1
             r_to   = PRC_RGN_lp2r(l_to,p_to)

             if( IO_L ) write(IO_FID_LOG,'(11(I6))') ipos,                 r_from, g_from, l_from, p_from, &
                                                           i_to  , j_to  , r_to  , g_to  , l_to  , p_to
          enddo
       enddo
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Recv_nmax_r2p(local)  = ', Recv_nmax_r2p
    if( IO_L ) write(IO_FID_LOG,*) '*** Send_nmax_r2p(local)  = ', Send_nmax_r2p
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(A)')             '|---------------------------------------'
    if( IO_L ) write(IO_FID_LOG,'(A)')             '|               size  prc_from    prc_to'
    if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))')    '| Copy_r2p', Copy_info_r2p(:)
    do irank = 1, Recv_nmax_r2p
       if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))') '| Recv_r2p', Recv_info_r2p(:,irank)
    enddo
    do irank = 1, Send_nmax_r2p
       if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))') '| Send_r2p', Send_info_r2p(:,irank)
    enddo

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Copy_list_r2p'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(11(A6))') 'number', &
                                               '|ifrom','|jfrom','|rfrom','|gfrom','|lfrom','|pfrom', &
                                                                 '|  rto','|  gto','|  lto','|  pto'
       do ipos = 1, Copy_info_r2p(I_size)
          g_from = Copy_list_r2p(I_grid_from,ipos)
          l_from = Copy_list_r2p(I_l_from,ipos)
          p_from = Copy_info_r2p(I_prc_from)
          i_from = mod(g_from-1,ADM_gall_1d) + 1
          j_from = (g_from-i_from) / ADM_gall_1d + 1
          r_from = PRC_RGN_lp2r(l_from,p_from)
          g_to   = Copy_list_r2p(I_grid_to,ipos)
          l_to   = Copy_list_r2p(I_l_to,ipos)
          p_to   = Copy_info_r2p(I_prc_to)
          r_to   = l_to

          if( IO_L ) write(IO_FID_LOG,'(11(I6))') ipos, i_from, j_from, r_from, g_from, l_from, p_from, &
                                                                        r_to  , g_to  , l_to  , p_to
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Recv_list_r2p'
       do irank = 1, Recv_nmax_r2p
          if( IO_L ) write(IO_FID_LOG,'(11(A6))') 'number', &
                                                  '|ifrom','|jfrom','|rfrom','|gfrom','|lfrom','|pfrom', &
                                                                    '|  rto','|  gto','|  lto','|  pto'
          do ipos = 1, Recv_info_r2p(I_size,irank)
             g_from = Recv_list_r2p(I_grid_from,ipos,irank)
             l_from = Recv_list_r2p(I_l_from,ipos,irank)
             p_from = Recv_info_r2p(I_prc_from,irank)
             i_from = mod(g_from-1,ADM_gall_1d) + 1
             j_from = (g_from-i_from) / ADM_gall_1d + 1
             r_from = PRC_RGN_lp2r(l_from,p_from)
             g_to   = Recv_list_r2p(I_grid_to,ipos,irank)
             l_to   = Recv_list_r2p(I_l_to,ipos,irank)
             p_to   = Recv_info_r2p(I_prc_to,irank)
             r_to   = l_to

             if( IO_L ) write(IO_FID_LOG,'(11(I6))') ipos, i_from, j_from, r_from, g_from, l_from, p_from, &
                                                                           r_to  , g_to  , l_to  , p_to
          enddo
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Send_list_r2p'
       do irank = 1, Send_nmax_r2p
          if( IO_L ) write(IO_FID_LOG,'(11(A6))') 'number', &
                                                  '|ifrom','|jfrom','|rfrom','|gfrom','|lfrom','|pfrom', &
                                                                    '|  rto','|  gto','|  lto','|  pto'
          do ipos = 1, Send_info_r2p(I_size,irank)
             g_from = Send_list_r2p(I_grid_from,ipos,irank)
             l_from = Send_list_r2p(I_l_from,ipos,irank)
             p_from = Send_info_r2p(I_prc_from,irank)
             i_from = mod(g_from-1,ADM_gall_1d) + 1
             j_from = (g_from-i_from) / ADM_gall_1d + 1
             r_from = PRC_RGN_lp2r(l_from,p_from)
             g_to   = Send_list_r2p(I_grid_to,ipos,irank)
             l_to   = Send_list_r2p(I_l_to,ipos,irank)
             p_to   = Send_info_r2p(I_prc_to,irank)
             r_to   = l_to

             if( IO_L ) write(IO_FID_LOG,'(11(I6))') ipos, i_from, j_from, r_from, g_from, l_from, p_from, &
                                                                           r_to  , g_to  , l_to  , p_to
          enddo
       enddo
    endif
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Send_size_p2r,r2p     = ', Send_size_nglobal_pl

    allocate( sendbuf_r2p_SP(Send_size_nglobal_pl*ADM_kall*COMM_varmax,Send_nmax_r2p) )
    allocate( recvbuf_r2p_SP(Send_size_nglobal_pl*ADM_kall*COMM_varmax,Recv_nmax_r2p) )
    allocate( sendbuf_p2r_SP(Send_size_nglobal_pl*ADM_kall*COMM_varmax,Send_nmax_p2r) )
    allocate( recvbuf_p2r_SP(Send_size_nglobal_pl*ADM_kall*COMM_varmax,Recv_nmax_p2r) )
    allocate( sendbuf_r2p_DP(Send_size_nglobal_pl*ADM_kall*COMM_varmax,Send_nmax_r2p) )
    allocate( recvbuf_r2p_DP(Send_size_nglobal_pl*ADM_kall*COMM_varmax,Recv_nmax_r2p) )
    allocate( sendbuf_p2r_DP(Send_size_nglobal_pl*ADM_kall*COMM_varmax,Send_nmax_p2r) )
    allocate( recvbuf_p2r_DP(Send_size_nglobal_pl*ADM_kall*COMM_varmax,Recv_nmax_p2r) )

    return
  end subroutine COMM_sortdest_pl

  !-----------------------------------------------------------------------------
  !> Setup data copy list of singular point
  subroutine COMM_sortdest_singular
    use scale_prc, only: &
       PRC_myrank
    use scale_prc_icoA, only: &
       PRC_RGN_local,    &
       PRC_RGN_vert_num, &
       PRC_RGN_lp2r,     &
       PRC_RGN_l2r,      &
       I_W,              &
       I_N,              &
       I_S
    implicit none

    integer  :: rgnid
    integer  :: i, j, i_rmt, j_rmt

    integer  :: l, ipos
    integer  :: i_from, j_from, r_from, g_from, l_from, p_from
    integer  :: i_to, j_to, r_to, g_to, l_to, p_to
    !---------------------------------------------------------------------------

    allocate( Singular_info(info_vindex) )
    Singular_info(:) = -1
    Singular_info(I_size) = 0

    allocate( Singular_list(list_vindex,4*ADM_lall) )
    Singular_list(:,:)   = -1

    do l = 1, PRC_RGN_local
       rgnid = PRC_RGN_l2r(l)

       if ( PRC_RGN_vert_num(I_W,rgnid) == 3 ) then
          Singular_info(I_size) = Singular_info(I_size) + 1
          ipos                  = Singular_info(I_size)

          i     = ADM_gmin
          j     = ADM_gmin - 1
          i_rmt = ADM_gmin - 1
          j_rmt = ADM_gmin - 1

          Singular_list(I_grid_from,ipos) = suf(i,j)
          Singular_list(I_l_from   ,ipos) = l
          Singular_list(I_grid_to  ,ipos) = suf(i_rmt,j_rmt)
          Singular_list(I_l_to     ,ipos) = l
       endif

       if ( PRC_RGN_vert_num(I_N,rgnid) /= 4 ) then
          Singular_info(I_size) = Singular_info(I_size) + 1
          ipos                  = Singular_info(I_size)

          i     = ADM_gmin
          j     = ADM_gmax + 1
          i_rmt = ADM_gmin - 1
          j_rmt = ADM_gmax + 1

          Singular_list(I_grid_from,ipos) = suf(i,j)
          Singular_list(I_l_from   ,ipos) = l
          Singular_list(I_grid_to  ,ipos) = suf(i_rmt,j_rmt)
          Singular_list(I_l_to     ,ipos) = l
       endif

       if ( PRC_RGN_vert_num(I_S,rgnid) /= 4 ) then
          Singular_info(I_size) = Singular_info(I_size) + 1
          ipos                  = Singular_info(I_size)

          i     = ADM_gmax + 1
          j     = ADM_gmin
          i_rmt = ADM_gmax + 1
          j_rmt = ADM_gmin - 1

          Singular_list(I_grid_from,ipos) = suf(i,j)
          Singular_list(I_l_from   ,ipos) = l
          Singular_list(I_grid_to  ,ipos) = suf(i_rmt,j_rmt)
          Singular_list(I_l_to     ,ipos) = l
       endif

    enddo ! loop l

    if ( Singular_info(I_size) > 0 ) then
       Singular_nmax = 1
       Singular_info(I_prc_from) = PRC_myrank
       Singular_info(I_prc_to  ) = PRC_myrank
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(A)')             '|---------------------------------------'
    if( IO_L ) write(IO_FID_LOG,'(A)')             '|               size  prc_from    prc_to'
    if( IO_L ) write(IO_FID_LOG,'(A10,3(I10))')    '| Singular', Singular_info(:)

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Singular_list'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(13(A6))') 'number', &
                                               '|ifrom','|jfrom','|rfrom','|gfrom','|lfrom','|pfrom', &
                                               '|  ito','|  jto','|  rto','|  gto','|  lto','|  pto'
       do ipos = 1, Singular_info(I_size)
          g_from = Singular_list(I_grid_from,ipos)
          l_from = Singular_list(I_l_from,ipos)
          p_from = Singular_info(I_prc_from)
          i_from = mod(g_from-1,ADM_gall_1d) + 1
          j_from = (g_from-i_from) / ADM_gall_1d + 1
          r_from = PRC_RGN_lp2r(l_from,p_from)
          g_to   = Singular_list(I_grid_to,ipos)
          l_to   = Singular_list(I_l_to,ipos)
          p_to   = Singular_info(I_prc_to)
          i_to   = mod(g_to-1,ADM_gall_1d) + 1
          j_to   = (g_to-i_to) / ADM_gall_1d + 1
          r_to   = PRC_RGN_lp2r(l_to,p_to)

          if( IO_L ) write(IO_FID_LOG,'(13(I6))') ipos, i_from, j_from, r_from, g_from, l_from, p_from, &
                                                        i_to  , j_to  , r_to  , g_to  , l_to  , p_to
       enddo
    endif

    return
  end subroutine COMM_sortdest_singular

  !-----------------------------------------------------------------------------
  !> Data transfer kernel
  subroutine COMM_data_transfer_SP( &
       var,   &
       var_pl )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_abort
    implicit none

    real(SP), intent(inout) :: var   (:,:,:,:)
    real(SP), intent(inout) :: var_pl(:,:,:,:)

    integer  :: shp(4), kmax, vmax
    integer  :: totalsize, rank, tag
    integer  :: irank, ipos, imax
    integer  :: ij_from, l_from, ij_to, l_to

    integer  :: k, v, ikv
    integer  :: ierr
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcopy(var) &
    !$acc pcopy(var_pl)

    if ( COMM_apply_barrier ) then
       call PROF_rapstart('COMM_barrier',2)
       call MPI_Barrier(PRC_LOCAL_COMM_WORLD,ierr)
       call PROF_rapend  ('COMM_barrier',2)
    endif

    !$acc wait

    call PROF_rapstart('COMM_data_transfer',2)

    shp  = shape(var)
    kmax = shp(2)
    vmax = shp(4)

    if ( kmax * vmax > ADM_kall * COMM_varmax ) then
       write(*,*) 'xxx [COMM_data_transfer] kmax * vmax exceeds ADM_kall * COMM_varmax, stop!'
       write(*,*) 'xxx kmax * vmax            = ', kmax * vmax
       write(*,*) 'xxx ADM_kall * COMM_varmax = ', ADM_kall * COMM_varmax
       call PRC_abort
    endif

    !---< start communication >---
    ! Theres no p2r & r2p communication without calling COMM_sortdest_pl.
    ! receive pole   => region
    ! receive region => pole
    ! receive region => region
    ! pack and send pole   => region
    ! pack and send region => pole
    ! pack and send region => region
    ! copy pole   => region
    ! copy region => pole
    ! copy region => region
    ! wait all
    ! unpack pole   => region
    ! unpack region => pole
    ! unpack region => region
    ! copy region halo => region halo (singular point)

    REQ_count = 0

    !--- receive r2r
    do irank = 1, Recv_nmax_r2r
       REQ_count = REQ_count + 1
       totalsize = Recv_info_r2r(I_size    ,irank) * kmax * vmax
       rank      = Recv_info_r2r(I_prc_from,irank)
       tag       = Recv_info_r2r(I_prc_from,irank)

       call MPI_IRECV( recvbuf_r2r_SP(1,irank), & ! [OUT]
                       totalsize,               & ! [IN]
                       MPI_REAL,                & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- receive p2r
    do irank = 1, Recv_nmax_p2r
       REQ_count = REQ_count + 1
       totalsize = Recv_info_p2r(I_size    ,irank) * kmax * vmax
       rank      = Recv_info_p2r(I_prc_from,irank)
       tag       = Recv_info_p2r(I_prc_from,irank) + 1000000

       call MPI_IRECV( recvbuf_p2r_SP(1,irank), & ! [OUT]
                       totalsize,               & ! [IN]
                       MPI_REAL,                & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- receive r2p
    do irank = 1, Recv_nmax_r2p
       REQ_count = REQ_count + 1
       totalsize = Recv_info_r2p(I_size    ,irank) * kmax * vmax
       rank      = Recv_info_r2p(I_prc_from,irank)
       tag       = Recv_info_r2p(I_prc_from,irank) + 2000000

       call MPI_IRECV( recvbuf_r2p_SP(1,irank), & ! [OUT]
                       totalsize,               & ! [IN]
                       MPI_REAL,                & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !$acc data &
    !$acc pcopyin(Send_list_r2r,Send_list_p2r,Send_list_r2p)
    !$acc wait

    !--- pack and send r2r
    do irank = 1, Send_nmax_r2r
       imax = Send_info_r2r(I_size,irank)

       !$acc kernels pcopy(sendbuf_r2r_SP)
       !$omp parallel do default(none),private(ipos,k,v,ikv,ij_from,l_from) &
       !$omp shared(irank,imax,kmax,vmax,Send_list_r2r,sendbuf_r2r_SP,var) &
       !$omp collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Send_list_r2r(I_grid_from,ipos,irank)
          l_from  = Send_list_r2r(I_l_from   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          sendbuf_r2r_SP(ikv,irank) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       REQ_count = REQ_count + 1
       totalsize = imax * kmax * vmax
       rank      = Send_info_r2r(I_prc_to  ,irank)
       tag       = Send_info_r2r(I_prc_from,irank)

       !$acc wait
       call MPI_ISEND( sendbuf_r2r_SP(1,irank), & ! [IN]
                       totalsize,               & ! [IN]
                       MPI_REAL,                & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- pack and send p2r
    do irank = 1, Send_nmax_p2r
       imax = Send_info_p2r(I_size,irank)

       !$acc kernels pcopy(sendbuf_p2r_SP)
       !$omp parallel do default(none),private(ipos,k,v,ikv,ij_from,l_from) &
       !$omp shared(irank,imax,kmax,vmax,Send_list_p2r,sendbuf_p2r_SP,var_pl) &
       !$omp collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Send_list_p2r(I_grid_from,ipos,irank)
          l_from  = Send_list_p2r(I_l_from   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          sendbuf_p2r_SP(ikv,irank) = var_pl(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       REQ_count = REQ_count + 1
       totalsize = imax * kmax * vmax
       rank      = Send_info_p2r(I_prc_to  ,irank)
       tag       = Send_info_p2r(I_prc_from,irank) + 1000000

       !$acc wait
       call MPI_ISEND( sendbuf_p2r_SP(1,irank), & ! [IN]
                       totalsize,               & ! [IN]
                       MPI_REAL,                & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- pack and send r2p
    do irank = 1, Send_nmax_r2p
       imax = Send_info_r2p(I_size,irank)

       !$acc kernels pcopy(sendbuf_r2p_SP)
       !$omp parallel do default(none),private(ipos,k,v,ikv,ij_from,l_from) &
       !$omp shared(irank,imax,kmax,vmax,Send_list_r2p,sendbuf_r2p_SP,var) &
       !$omp collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Send_list_r2p(I_grid_from,ipos,irank)
          l_from  = Send_list_r2p(I_l_from   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          sendbuf_r2p_SP(ikv,irank) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       REQ_count = REQ_count + 1
       totalsize = imax * kmax * vmax
       rank      = Send_info_r2p(I_prc_to  ,irank)
       tag       = Send_info_r2p(I_prc_from,irank) + 2000000

       !$acc wait
       call MPI_ISEND( sendbuf_r2p_SP(1,irank), & ! [IN]
                       totalsize,               & ! [IN]
                       MPI_REAL,                & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !$acc end data
    !$acc data &
    !$acc pcopyin(Copy_list_r2r,Copy_list_p2r,Copy_list_r2p)
    !$acc wait

    !$omp parallel default(none),private(ipos,k,v,imax,irank,ij_from,l_from,ij_to,l_to) &
    !$omp shared(kmax,vmax,var,var_pl,                      &
    !$omp        Copy_nmax_p2r,Copy_info_p2r,Copy_list_p2r, &
    !$omp        Copy_nmax_r2p,Copy_info_r2p,Copy_list_r2p, &
    !$omp        Copy_nmax_r2r,Copy_info_r2r,Copy_list_r2r  )

    !--- copy r2r
    do irank = 1, Copy_nmax_r2r
       imax = Copy_info_r2r(I_size)

       !$acc kernels
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Copy_list_r2r(I_grid_from,ipos)
          l_from  = Copy_list_r2r(I_l_from   ,ipos)
          ij_to   = Copy_list_r2r(I_grid_to  ,ipos)
          l_to    = Copy_list_r2r(I_l_to     ,ipos)

          var(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !--- copy p2r
    do irank = 1, Copy_nmax_p2r
       imax = Copy_info_p2r(I_size)

       !$acc kernels
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Copy_list_p2r(I_grid_from,ipos)
          l_from  = Copy_list_p2r(I_l_from   ,ipos)
          ij_to   = Copy_list_p2r(I_grid_to  ,ipos)
          l_to    = Copy_list_p2r(I_l_to     ,ipos)

          var(ij_to,k,l_to,v) = var_pl(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !--- copy r2p
    do irank = 1, Copy_nmax_r2p
       imax = Copy_info_r2p(I_size)

       !$acc kernels
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Copy_list_r2p(I_grid_from,ipos)
          l_from  = Copy_list_r2p(I_l_from   ,ipos)
          ij_to   = Copy_list_r2p(I_grid_to  ,ipos)
          l_to    = Copy_list_r2p(I_l_to     ,ipos)

          var_pl(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !$omp end parallel

    !--- wait all
    if ( REQ_count > 0 ) then
       call MPI_WAITALL( REQ_count,           & ! [IN]
                         REQ_list(:),         & ! [IN]
                         MPI_STATUSES_IGNORE, & ! [OUT]
                         ierr                 ) ! [OUT]
    endif

    !$acc end data
    !$acc data &
    !$acc pcopyin(Recv_list_r2r,Recv_list_p2r,Recv_list_r2p)
    !$acc wait

    !$omp parallel default(none),private(ipos,k,v,imax,irank,ikv,ij_from,l_from,ij_to,l_to) &
    !$omp shared(kmax,vmax,var,var_pl,                                     &
    !$omp        Recv_nmax_p2r,Recv_info_p2r,Recv_list_p2r,recvbuf_p2r_SP, &
    !$omp        Recv_nmax_r2p,Recv_info_r2p,Recv_list_r2p,recvbuf_r2p_SP, &
    !$omp        Recv_nmax_r2r,Recv_info_r2r,Recv_list_r2r,recvbuf_r2r_SP, &
    !$omp        Singular_nmax,Singular_info,Singular_list                 )

    !--- unpack r2r
    do irank = 1, Recv_nmax_r2r
       imax = Recv_info_r2r(I_size,irank)

       !$acc kernels pcopyin(recvbuf_r2r_SP)
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_to = Recv_list_r2r(I_grid_to,ipos,irank)
          l_to  = Recv_list_r2r(I_l_to   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          var(ij_to,k,l_to,v) = recvbuf_r2r_SP(ikv,irank)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !--- unpack p2r
    do irank = 1, Recv_nmax_p2r
       imax = Recv_info_p2r(I_size,irank)

       !$acc kernels pcopyin(recvbuf_p2r_SP)
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_to = Recv_list_p2r(I_grid_to,ipos,irank)
          l_to  = Recv_list_p2r(I_l_to   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          var(ij_to,k,l_to,v) = recvbuf_p2r_SP(ikv,irank)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !--- unpack r2p
    do irank = 1, Recv_nmax_r2p
       imax = Recv_info_r2p(I_size,irank)

       !$acc kernels pcopyin(recvbuf_r2p_SP)
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_to = Recv_list_r2p(I_grid_to,ipos,irank)
          l_to  = Recv_list_r2p(I_l_to   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          var_pl(ij_to,k,l_to,v) = recvbuf_r2p_SP(ikv,irank)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !$acc end data
    !$acc wait

    !--- singular point (halo to halo)
    do irank = 1, Singular_nmax
       imax = Singular_info(I_size)

       !$acc kernels pcopyin(Singular_list)
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Singular_list(I_grid_from,ipos)
          l_from  = Singular_list(I_l_from   ,ipos)
          ij_to   = Singular_list(I_grid_to  ,ipos)
          l_to    = Singular_list(I_l_to     ,ipos)

          var(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !$omp end parallel

    call PROF_rapend('COMM_data_transfer',2)

    !$acc end data
    !$acc wait

    return
  end subroutine COMM_data_transfer_SP

  !-----------------------------------------------------------------------------
  !> Data transfer kernel
  subroutine COMM_data_transfer_DP( &
       var,   &
       var_pl )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_abort
    implicit none

    real(DP), intent(inout) :: var   (:,:,:,:)
    real(DP), intent(inout) :: var_pl(:,:,:,:)

    integer  :: shp(4), kmax, vmax
    integer  :: totalsize, rank, tag
    integer  :: irank, ipos, imax
    integer  :: ij_from, l_from, ij_to, l_to

    integer  :: k, v, ikv
    integer  :: ierr
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcopy(var) &
    !$acc pcopy(var_pl)

    if ( COMM_apply_barrier ) then
       call PROF_rapstart('COMM_barrier',2)
       call MPI_Barrier(PRC_LOCAL_COMM_WORLD,ierr)
       call PROF_rapend  ('COMM_barrier',2)
    endif

    !$acc wait

    call PROF_rapstart('COMM_data_transfer',2)

    shp  = shape(var)
    kmax = shp(2)
    vmax = shp(4)

    if ( kmax * vmax > ADM_kall * COMM_varmax ) then
       write(*,*) 'xxx [COMM_data_transfer] kmax * vmax exceeds ADM_kall * COMM_varmax, stop!'
       write(*,*) 'xxx kmax * vmax            = ', kmax * vmax
       write(*,*) 'xxx ADM_kall * COMM_varmax = ', ADM_kall * COMM_varmax
       call PRC_abort
    endif

    !---< start communication >---
    ! Theres no p2r & r2p communication without calling COMM_sortdest_pl.
    ! receive pole   => region
    ! receive region => pole
    ! receive region => region
    ! pack and send pole   => region
    ! pack and send region => pole
    ! pack and send region => region
    ! copy pole   => region
    ! copy region => pole
    ! copy region => region
    ! wait all
    ! unpack pole   => region
    ! unpack region => pole
    ! unpack region => region
    ! copy region halo => region halo (singular point)

    REQ_count = 0

    !--- receive r2r
    do irank = 1, Recv_nmax_r2r
       REQ_count = REQ_count + 1
       totalsize = Recv_info_r2r(I_size    ,irank) * kmax * vmax
       rank      = Recv_info_r2r(I_prc_from,irank)
       tag       = Recv_info_r2r(I_prc_from,irank)

       call MPI_IRECV( recvbuf_r2r_DP(1,irank), & ! [OUT]
                       totalsize,               & ! [IN]
                       MPI_DOUBLE_PRECISION,    & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- receive p2r
    do irank = 1, Recv_nmax_p2r
       REQ_count = REQ_count + 1
       totalsize = Recv_info_p2r(I_size    ,irank) * kmax * vmax
       rank      = Recv_info_p2r(I_prc_from,irank)
       tag       = Recv_info_p2r(I_prc_from,irank) + 1000000

       call MPI_IRECV( recvbuf_p2r_DP(1,irank), & ! [OUT]
                       totalsize,               & ! [IN]
                       MPI_DOUBLE_PRECISION,    & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- receive r2p
    do irank = 1, Recv_nmax_r2p
       REQ_count = REQ_count + 1
       totalsize = Recv_info_r2p(I_size    ,irank) * kmax * vmax
       rank      = Recv_info_r2p(I_prc_from,irank)
       tag       = Recv_info_r2p(I_prc_from,irank) + 2000000

       call MPI_IRECV( recvbuf_r2p_DP(1,irank), & ! [OUT]
                       totalsize,               & ! [IN]
                       MPI_DOUBLE_PRECISION,    & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !$acc data &
    !$acc pcopyin(Send_list_r2r,Send_list_p2r,Send_list_r2p)
    !$acc wait

    !--- pack and send r2r
    do irank = 1, Send_nmax_r2r
       imax = Send_info_r2r(I_size,irank)

       !$acc kernels pcopy(sendbuf_r2r_DP)
       !$omp parallel do default(none),private(ipos,k,v,ikv,ij_from,l_from) &
       !$omp shared(irank,imax,kmax,vmax,Send_list_r2r,sendbuf_r2r_DP,var) &
       !$omp collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Send_list_r2r(I_grid_from,ipos,irank)
          l_from  = Send_list_r2r(I_l_from   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          sendbuf_r2r_DP(ikv,irank) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       REQ_count = REQ_count + 1
       totalsize = imax * kmax * vmax
       rank      = Send_info_r2r(I_prc_to  ,irank)
       tag       = Send_info_r2r(I_prc_from,irank)

       !$acc wait
       call MPI_ISEND( sendbuf_r2r_DP(1,irank), & ! [IN]
                       totalsize,               & ! [IN]
                       MPI_DOUBLE_PRECISION,    & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- pack and send p2r
    do irank = 1, Send_nmax_p2r
       imax = Send_info_p2r(I_size,irank)

       !$acc kernels pcopy(sendbuf_p2r_DP)
       !$omp parallel do default(none),private(ipos,k,v,ikv,ij_from,l_from) &
       !$omp shared(irank,imax,kmax,vmax,Send_list_p2r,sendbuf_p2r_DP,var_pl) &
       !$omp collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Send_list_p2r(I_grid_from,ipos,irank)
          l_from  = Send_list_p2r(I_l_from   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          sendbuf_p2r_DP(ikv,irank) = var_pl(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       REQ_count = REQ_count + 1
       totalsize = imax * kmax * vmax
       rank      = Send_info_p2r(I_prc_to  ,irank)
       tag       = Send_info_p2r(I_prc_from,irank) + 1000000

       !$acc wait
       call MPI_ISEND( sendbuf_p2r_DP(1,irank), & ! [IN]
                       totalsize,               & ! [IN]
                       MPI_DOUBLE_PRECISION,    & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- pack and send r2p
    do irank = 1, Send_nmax_r2p
       imax = Send_info_r2p(I_size,irank)

       !$acc kernels pcopy(sendbuf_r2p_DP)
       !$omp parallel do default(none),private(ipos,k,v,ikv,ij_from,l_from) &
       !$omp shared(irank,imax,kmax,vmax,Send_list_r2p,sendbuf_r2p_DP,var) &
       !$omp collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Send_list_r2p(I_grid_from,ipos,irank)
          l_from  = Send_list_r2p(I_l_from   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          sendbuf_r2p_DP(ikv,irank) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       REQ_count = REQ_count + 1
       totalsize = imax * kmax * vmax
       rank      = Send_info_r2p(I_prc_to  ,irank)
       tag       = Send_info_r2p(I_prc_from,irank) + 2000000

       !$acc wait
       call MPI_ISEND( sendbuf_r2p_DP(1,irank), & ! [IN]
                       totalsize,               & ! [IN]
                       MPI_DOUBLE_PRECISION,    & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !$acc end data
    !$acc data &
    !$acc pcopyin(Copy_list_r2r,Copy_list_p2r,Copy_list_r2p)
    !$acc wait

    !$omp parallel default(none),private(ipos,k,v,imax,irank,ij_from,l_from,ij_to,l_to) &
    !$omp shared(kmax,vmax,var,var_pl,                      &
    !$omp        Copy_nmax_p2r,Copy_info_p2r,Copy_list_p2r, &
    !$omp        Copy_nmax_r2p,Copy_info_r2p,Copy_list_r2p, &
    !$omp        Copy_nmax_r2r,Copy_info_r2r,Copy_list_r2r  )

    !--- copy r2r
    do irank = 1, Copy_nmax_r2r
       imax = Copy_info_r2r(I_size)

       !$acc kernels
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Copy_list_r2r(I_grid_from,ipos)
          l_from  = Copy_list_r2r(I_l_from   ,ipos)
          ij_to   = Copy_list_r2r(I_grid_to  ,ipos)
          l_to    = Copy_list_r2r(I_l_to     ,ipos)

          var(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !--- copy p2r
    do irank = 1, Copy_nmax_p2r
       imax = Copy_info_p2r(I_size)

       !$acc kernels
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Copy_list_p2r(I_grid_from,ipos)
          l_from  = Copy_list_p2r(I_l_from   ,ipos)
          ij_to   = Copy_list_p2r(I_grid_to  ,ipos)
          l_to    = Copy_list_p2r(I_l_to     ,ipos)

          var(ij_to,k,l_to,v) = var_pl(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !--- copy r2p
    do irank = 1, Copy_nmax_r2p
       imax = Copy_info_r2p(I_size)

       !$acc kernels
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Copy_list_r2p(I_grid_from,ipos)
          l_from  = Copy_list_r2p(I_l_from   ,ipos)
          ij_to   = Copy_list_r2p(I_grid_to  ,ipos)
          l_to    = Copy_list_r2p(I_l_to     ,ipos)

          var_pl(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !$omp end parallel

    !--- wait all
    if ( REQ_count > 0 ) then
       call MPI_WAITALL( REQ_count,           & ! [IN]
                         REQ_list(:),         & ! [IN]
                         MPI_STATUSES_IGNORE, & ! [OUT]
                         ierr                 ) ! [OUT]
    endif

    !$acc end data
    !$acc data &
    !$acc pcopyin(Recv_list_r2r,Recv_list_p2r,Recv_list_r2p)
    !$acc wait

    !$omp parallel default(none),private(ipos,k,v,imax,irank,ikv,ij_from,l_from,ij_to,l_to) &
    !$omp shared(kmax,vmax,var,var_pl,                                     &
    !$omp        Recv_nmax_p2r,Recv_info_p2r,Recv_list_p2r,recvbuf_p2r_DP, &
    !$omp        Recv_nmax_r2p,Recv_info_r2p,Recv_list_r2p,recvbuf_r2p_DP, &
    !$omp        Recv_nmax_r2r,Recv_info_r2r,Recv_list_r2r,recvbuf_r2r_DP, &
    !$omp        Singular_nmax,Singular_info,Singular_list                 )

    !--- unpack r2r
    do irank = 1, Recv_nmax_r2r
       imax = Recv_info_r2r(I_size,irank)

       !$acc kernels pcopyin(recvbuf_r2r_DP)
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_to = Recv_list_r2r(I_grid_to,ipos,irank)
          l_to  = Recv_list_r2r(I_l_to   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          var(ij_to,k,l_to,v) = recvbuf_r2r_DP(ikv,irank)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !--- unpack p2r
    do irank = 1, Recv_nmax_p2r
       imax = Recv_info_p2r(I_size,irank)

       !$acc kernels pcopyin(recvbuf_p2r_DP)
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_to = Recv_list_p2r(I_grid_to,ipos,irank)
          l_to  = Recv_list_p2r(I_l_to   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          var(ij_to,k,l_to,v) = recvbuf_p2r_DP(ikv,irank)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !--- unpack r2p
    do irank = 1, Recv_nmax_r2p
       imax = Recv_info_r2p(I_size,irank)

       !$acc kernels pcopyin(recvbuf_r2p_DP)
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_to = Recv_list_r2p(I_grid_to,ipos,irank)
          l_to  = Recv_list_r2p(I_l_to   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          var_pl(ij_to,k,l_to,v) = recvbuf_r2p_DP(ikv,irank)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !$acc end data
    !$acc wait

    !--- singular point (halo to halo)
    do irank = 1, Singular_nmax
       imax = Singular_info(I_size)

       !$acc kernels pcopyin(Singular_list)
       !$omp do collapse(2)
       !$acc loop independent
       do v    = 1, vmax
       !$acc loop independent
       do k    = 1, kmax
       !$acc loop independent
       do ipos = 1, imax
          ij_from = Singular_list(I_grid_from,ipos)
          l_from  = Singular_list(I_l_from   ,ipos)
          ij_to   = Singular_list(I_grid_to  ,ipos)
          l_to    = Singular_list(I_l_to     ,ipos)

          var(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$acc end kernels
    enddo

    !$omp end parallel

    call PROF_rapend('COMM_data_transfer',2)

    !$acc end data
    !$acc wait

    return
  end subroutine COMM_data_transfer_DP

  !-----------------------------------------------------------------------------
  !> Data transfer kernel
  subroutine COMM_data_transfer_nopl( &
       var )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_abort
    implicit none

    real(DP), intent(inout) :: var   (:,:,:,:)

    integer  :: shp(4), kmax, vmax
    integer  :: totalsize, rank, tag
    integer  :: irank, ipos, imax
    integer  :: ij_from, l_from, ij_to, l_to

    integer  :: k, v, ikv
    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_apply_barrier ) then
       call PROF_rapstart('COMM_barrier',2)
       call MPI_Barrier(PRC_LOCAL_COMM_WORLD,ierr)
       call PROF_rapend  ('COMM_barrier',2)
    endif

    call PROF_rapstart('COMM_data_transfer',2)

    shp  = shape(var)
    kmax = shp(2)
    vmax = shp(4)

    if ( kmax * vmax > ADM_kall * COMM_varmax ) then
       write(*,*) 'xxx [COMM_data_transfer_nopl] kmax * vmax exceeds ADM_kall * COMM_varmax, stop!'
       write(*,*) 'xxx kmax * vmax            = ', kmax * vmax
       write(*,*) 'xxx ADM_kall * COMM_varmax = ', ADM_kall * COMM_varmax
       call PRC_abort
    endif

    !---< start communication >---
    ! receive region => region
    ! pack and send region => region
    ! copy region => region
    ! wait all
    ! unpack region => region
    ! copy region halo => region halo (singular point)

    REQ_count = 0

    !--- receive r2r
    do irank = 1, Recv_nmax_r2r
       REQ_count = REQ_count + 1
       totalsize = Recv_info_r2r(I_size    ,irank) * kmax * vmax
       rank      = Recv_info_r2r(I_prc_from,irank)
       tag       = Recv_info_r2r(I_prc_from,irank)

       call MPI_IRECV( recvbuf_r2r_DP(1,irank), & ! [OUT]
                       totalsize,               & ! [IN]
                       MPI_DOUBLE_PRECISION,    & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- pack and send r2r
    do irank = 1, Send_nmax_r2r
       imax = Send_info_r2r(I_size,irank)

       !$omp parallel do default(none),private(ipos,k,v,ikv,ij_from,l_from) &
       !$omp shared(irank,imax,kmax,vmax,Send_list_r2r,sendbuf_r2r_DP,var)
       do v    = 1, vmax
       do k    = 1, kmax
       do ipos = 1, imax
          ij_from = Send_list_r2r(I_grid_from,ipos,irank)
          l_from  = Send_list_r2r(I_l_from   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          sendbuf_r2r_DP(ikv,irank) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
       !$omp end parallel do

       REQ_count = REQ_count + 1
       totalsize = imax * kmax * vmax
       rank      = Send_info_r2r(I_prc_to  ,irank)
       tag       = Send_info_r2r(I_prc_from,irank)

       call MPI_ISEND( sendbuf_r2r_DP(1,irank), & ! [IN]
                       totalsize,               & ! [IN]
                       MPI_DOUBLE_PRECISION,    & ! [IN]
                       rank,                    & ! [IN]
                       tag,                     & ! [IN]
                       PRC_LOCAL_COMM_WORLD,    & ! [IN]
                       REQ_list(REQ_count),     & ! [OUT]
                       ierr                     ) ! [OUT]
    enddo

    !--- copy r2r
    !$omp parallel do default(none),private(ipos,k,v,imax,irank,ij_from,l_from,ij_to,l_to) &
    !$omp shared(kmax,vmax,Copy_nmax_r2r,Copy_info_r2r,Copy_list_r2r,var)
    do irank = 1, Copy_nmax_r2r
       imax = Copy_info_r2r(I_size)

       do v    = 1, vmax
       do k    = 1, kmax
       do ipos = 1, imax
          ij_from = Copy_list_r2r(I_grid_from,ipos)
          l_from  = Copy_list_r2r(I_l_from   ,ipos)
          ij_to   = Copy_list_r2r(I_grid_to  ,ipos)
          l_to    = Copy_list_r2r(I_l_to     ,ipos)

          var(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
    enddo
    !$omp end parallel do

    !--- wait all
    if ( REQ_count > 0 ) then
       call MPI_WAITALL( REQ_count,           & ! [IN]
                         REQ_list(:),         & ! [IN]
                         MPI_STATUSES_IGNORE, & ! [OUT]
                         ierr                 ) ! [OUT]
    endif

    !--- unpack r2r
    !$omp parallel do default(none),private(ipos,k,v,imax,irank,ij_to,l_to,ikv) &
    !$omp shared(kmax,vmax,Recv_nmax_r2r,Recv_info_r2r,Recv_list_r2r,recvbuf_r2r_DP,var)
    do irank = 1, Recv_nmax_r2r
       imax = Recv_info_r2r(I_size,irank)

       do v    = 1, vmax
       do k    = 1, kmax
       do ipos = 1, imax
          ij_to = Recv_list_r2r(I_grid_to,ipos,irank)
          l_to  = Recv_list_r2r(I_l_to   ,ipos,irank)

          ikv = (v-1) * imax * kmax &
              + (k-1) * imax        &
              + ipos

          var(ij_to,k,l_to,v) = recvbuf_r2r_DP(ikv,irank)
       enddo
       enddo
       enddo
    enddo
    !$omp end parallel do

    !--- singular point (halo to halo)
    do irank = 1, Singular_nmax
       imax = Singular_info(I_size)

       do v    = 1, vmax
       do k    = 1, kmax
       do ipos = 1, imax
          ij_from = Singular_list(I_grid_from,ipos)
          l_from  = Singular_list(I_l_from   ,ipos)
          ij_to   = Singular_list(I_grid_to  ,ipos)
          l_to    = Singular_list(I_l_to     ,ipos)

          var(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
       enddo
       enddo
       enddo
    enddo

    call PROF_rapend('COMM_data_transfer',2)

    return
  end subroutine COMM_data_transfer_nopl

  !-----------------------------------------------------------------------------
  !> Data transfer with region halo => pole center
  subroutine COMM_var_SP( &
       var,    &
       var_pl, &
       kmax,   &
       vmax    )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD
    use scale_prc_icoA, only: &
       PRC_RGN_total_pl, &
       PRC_RGN_rgn4pl,   &
       PRC_RGN_lp2r,     &
       I_NPL,            &
       I_SPL
    implicit none

    integer,  intent(in)    :: kmax
    integer,  intent(in)    :: vmax
    real(SP), intent(inout) :: var   (ADM_gall,   kmax,ADM_lall,   vmax)
    real(SP), intent(inout) :: var_pl(ADM_gall_pl,kmax,ADM_lall_pl,vmax)

    real(SP) :: sendbuf_h2p_SP(kmax*vmax,PRC_RGN_total_pl)
    real(SP) :: recvbuf_h2p_SP(kmax*vmax,PRC_RGN_total_pl)

    integer  :: totalsize, rank, tag
    integer  :: irank, ipos
    integer  :: ij_from, l_from, ij_to, l_to
    integer  :: r_from, r_to

    integer  :: k, v, kk
    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_apply_barrier ) then
       call PROF_rapstart('COMM_barrier',2)
       call MPI_Barrier(PRC_LOCAL_COMM_WORLD,ierr)
       call PROF_rapend  ('COMM_barrier',2)
    endif

    call PROF_rapstart('COMM_var',2)

    if ( comm_pl ) then

    REQ_count = 0

    !--- receive p2r-reverse
    do irank = 1, Send_nmax_p2r
       do ipos = 1, Send_info_p2r(I_size,irank)
          l_from = Send_list_p2r(I_l_to   ,ipos,irank)
          r_from = PRC_RGN_lp2r(l_from,Send_info_p2r(I_prc_to,irank))

          if ( r_from == PRC_RGN_rgn4pl(I_NPL) ) then
             REQ_count = REQ_count + 1
             totalsize = kmax * vmax
             rank      = Send_info_p2r(I_prc_to  ,irank)
             tag       = Send_info_p2r(I_prc_from,irank) + 1000000

             call MPI_IRECV( recvbuf_h2p_SP(1,I_NPL), & ! [OUT]
                             totalsize,               & ! [IN]
                             MPI_REAL,                & ! [IN]
                             rank,                    & ! [IN]
                             tag,                     & ! [IN]
                             PRC_LOCAL_COMM_WORLD,    & ! [IN]
                             REQ_list(REQ_count),     & ! [OUT]
                             ierr                     ) ! [OUT]
          endif

          if ( r_from == PRC_RGN_rgn4pl(I_SPL) ) then
             REQ_count = REQ_count + 1
             totalsize = kmax * vmax
             rank      = Send_info_p2r(I_prc_to  ,irank)
             tag       = Send_info_p2r(I_prc_from,irank) + 2000000

             call MPI_IRECV( recvbuf_h2p_SP(1,I_SPL), & ! [OUT]
                             totalsize,               & ! [IN]
                             MPI_REAL,                & ! [IN]
                             rank,                    & ! [IN]
                             tag,                     & ! [IN]
                             PRC_LOCAL_COMM_WORLD,    & ! [IN]
                             REQ_list(REQ_count),     & ! [OUT]
                             ierr                     ) ! [OUT]
          endif
       enddo
    enddo

    !--- pack and send p2r-reverse
    do irank = 1, Recv_nmax_p2r
       do ipos = 1, Recv_info_p2r(I_size,irank)
          ij_from = Recv_list_p2r(I_grid_to,ipos,irank)
          l_from  = Recv_list_p2r(I_l_to   ,ipos,irank)
          r_from  = PRC_RGN_lp2r(l_from,Recv_info_p2r(I_prc_to,irank))

          if ( r_from == PRC_RGN_rgn4pl(I_NPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                kk = (v-1) * kmax + k
                sendbuf_h2p_SP(kk,I_NPL) = var(ij_from,k,l_from,v)
             enddo
             enddo

             REQ_count = REQ_count + 1
             totalsize = kmax * vmax
             rank      = Recv_info_p2r(I_prc_from,irank)
             tag       = Recv_info_p2r(I_prc_from,irank) + 1000000

             call MPI_ISEND( sendbuf_h2p_SP(1,I_NPL), & ! [IN]
                             totalsize,               & ! [IN]
                             MPI_REAL,                & ! [IN]
                             rank,                    & ! [IN]
                             tag,                     & ! [IN]
                             PRC_LOCAL_COMM_WORLD,    & ! [IN]
                             REQ_list(REQ_count),     & ! [OUT]
                             ierr                     ) ! [OUT]
          endif

          if ( r_from == PRC_RGN_rgn4pl(I_SPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                kk = (v-1) * kmax + k
                sendbuf_h2p_SP(kk,I_SPL) = var(ij_from,k,l_from,v)
             enddo
             enddo

             REQ_count = REQ_count + 1
             totalsize = kmax * vmax
             rank      = Recv_info_p2r(I_prc_from,irank)
             tag       = Recv_info_p2r(I_prc_from,irank) + 2000000

             call MPI_ISEND( sendbuf_h2p_SP(1,I_SPL), & ! [IN]
                             totalsize,               & ! [IN]
                             MPI_REAL,                & ! [IN]
                             rank,                    & ! [IN]
                             tag,                     & ! [IN]
                             PRC_LOCAL_COMM_WORLD,    & ! [IN]
                             REQ_list(REQ_count),     & ! [OUT]
                             ierr                     ) ! [OUT]
          endif
       enddo
    enddo

    !--- copy p2r-reverse
    do irank = 1, Copy_nmax_p2r
       do ipos = 1, Copy_info_p2r(I_size)
          ij_from = Copy_list_p2r(I_grid_to  ,ipos)
          l_from  = Copy_list_p2r(I_l_to     ,ipos)
          r_from  = PRC_RGN_lp2r(l_from,Copy_info_p2r(I_prc_to))
          ij_to   = Copy_list_p2r(I_grid_from,ipos)
          l_to    = Copy_list_p2r(I_l_from   ,ipos)
          r_to    = PRC_RGN_lp2r(l_to,Copy_info_p2r(I_prc_from))

          if (      r_from == PRC_RGN_rgn4pl(I_NPL) &
               .OR. r_from == PRC_RGN_rgn4pl(I_SPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                var_pl(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
             enddo
             enddo
          endif
       enddo
    enddo

    !--- wait all
    if ( REQ_count > 0 ) then
       call MPI_WAITALL( REQ_count,           & ! [IN]
                         REQ_list(:),         & ! [IN]
                         MPI_STATUSES_IGNORE, & ! [OUT]
                         ierr                 ) ! [OUT]
    endif

    !--- unpack p2r-reverse
    do irank = 1, Send_nmax_p2r
       do ipos = 1, Send_info_p2r(I_size,irank)
          l_from = Send_list_p2r(I_l_to     ,ipos,irank)
          r_from = PRC_RGN_lp2r(l_from,Send_info_p2r(I_prc_to,irank))

          ij_to  = Send_list_p2r(I_grid_from,ipos,irank)
          l_to   = Send_list_p2r(I_l_from   ,ipos,irank)

          if ( r_from == PRC_RGN_rgn4pl(I_NPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                kk = (v-1) * kmax + k
                var_pl(ij_to,k,l_to,v) = recvbuf_h2p_SP(kk,I_NPL)
             enddo
             enddo
          endif

          if ( r_from == PRC_RGN_rgn4pl(I_SPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                kk = (v-1) * kmax + k
                var_pl(ij_to,k,l_to,v) = recvbuf_h2p_SP(kk,I_SPL)
             enddo
             enddo
          endif
       enddo
    enddo

    endif

    call COMM_data_transfer_SP(var,var_pl)

    call PROF_rapend('COMM_var',2)

    return
  end subroutine COMM_var_SP

  !-----------------------------------------------------------------------------
  !> Data transfer with region halo => pole center
  subroutine COMM_var_DP( &
       var,    &
       var_pl, &
       kmax,   &
       vmax    )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD
    use scale_prc_icoA, only: &
       PRC_RGN_total_pl, &
       PRC_RGN_rgn4pl,   &
       PRC_RGN_lp2r,     &
       I_NPL,            &
       I_SPL
    implicit none

    integer,  intent(in)    :: kmax
    integer,  intent(in)    :: vmax
    real(DP), intent(inout) :: var   (ADM_gall,   kmax,ADM_lall,   vmax)
    real(DP), intent(inout) :: var_pl(ADM_gall_pl,kmax,ADM_lall_pl,vmax)

    real(DP) :: sendbuf_h2p_DP(kmax*vmax,PRC_RGN_total_pl)
    real(DP) :: recvbuf_h2p_DP(kmax*vmax,PRC_RGN_total_pl)

    integer  :: totalsize, rank, tag
    integer  :: irank, ipos
    integer  :: ij_from, l_from, ij_to, l_to
    integer  :: r_from, r_to

    integer  :: k, v, kk
    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_apply_barrier ) then
       call PROF_rapstart('COMM_barrier',2)
       call MPI_Barrier(PRC_LOCAL_COMM_WORLD,ierr)
       call PROF_rapend  ('COMM_barrier',2)
    endif

    call PROF_rapstart('COMM_var',2)

    if ( comm_pl ) then

    REQ_count = 0

    !--- receive p2r-reverse
    do irank = 1, Send_nmax_p2r
       do ipos = 1, Send_info_p2r(I_size,irank)
          l_from = Send_list_p2r(I_l_to   ,ipos,irank)
          r_from = PRC_RGN_lp2r(l_from,Send_info_p2r(I_prc_to,irank))

          if ( r_from == PRC_RGN_rgn4pl(I_NPL) ) then
             REQ_count = REQ_count + 1
             totalsize = kmax * vmax
             rank      = Send_info_p2r(I_prc_to  ,irank)
             tag       = Send_info_p2r(I_prc_from,irank) + 1000000

             call MPI_IRECV( recvbuf_h2p_DP(1,I_NPL), & ! [OUT]
                             totalsize,               & ! [IN]
                             MPI_DOUBLE_PRECISION,    & ! [IN]
                             rank,                    & ! [IN]
                             tag,                     & ! [IN]
                             PRC_LOCAL_COMM_WORLD,    & ! [IN]
                             REQ_list(REQ_count),     & ! [OUT]
                             ierr                     ) ! [OUT]
          endif

          if ( r_from == PRC_RGN_rgn4pl(I_SPL) ) then
             REQ_count = REQ_count + 1
             totalsize = kmax * vmax
             rank      = Send_info_p2r(I_prc_to  ,irank)
             tag       = Send_info_p2r(I_prc_from,irank) + 2000000

             call MPI_IRECV( recvbuf_h2p_DP(1,I_SPL), & ! [OUT]
                             totalsize,               & ! [IN]
                             MPI_DOUBLE_PRECISION,    & ! [IN]
                             rank,                    & ! [IN]
                             tag,                     & ! [IN]
                             PRC_LOCAL_COMM_WORLD,    & ! [IN]
                             REQ_list(REQ_count),     & ! [OUT]
                             ierr                     ) ! [OUT]
          endif
       enddo
    enddo

    !--- pack and send p2r-reverse
    do irank = 1, Recv_nmax_p2r
       do ipos = 1, Recv_info_p2r(I_size,irank)
          ij_from = Recv_list_p2r(I_grid_to,ipos,irank)
          l_from  = Recv_list_p2r(I_l_to   ,ipos,irank)
          r_from  = PRC_RGN_lp2r(l_from,Recv_info_p2r(I_prc_to,irank))

          if ( r_from == PRC_RGN_rgn4pl(I_NPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                kk = (v-1) * kmax + k
                sendbuf_h2p_DP(kk,I_NPL) = var(ij_from,k,l_from,v)
             enddo
             enddo

             REQ_count = REQ_count + 1
             totalsize = kmax * vmax
             rank      = Recv_info_p2r(I_prc_from,irank)
             tag       = Recv_info_p2r(I_prc_from,irank) + 1000000

             call MPI_ISEND( sendbuf_h2p_DP(1,I_NPL), & ! [IN]
                             totalsize,               & ! [IN]
                             MPI_DOUBLE_PRECISION,    & ! [IN]
                             rank,                    & ! [IN]
                             tag,                     & ! [IN]
                             PRC_LOCAL_COMM_WORLD,    & ! [IN]
                             REQ_list(REQ_count),     & ! [OUT]
                             ierr                     ) ! [OUT]
          endif

          if ( r_from == PRC_RGN_rgn4pl(I_SPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                kk = (v-1) * kmax + k
                sendbuf_h2p_DP(kk,I_SPL) = var(ij_from,k,l_from,v)
             enddo
             enddo

             REQ_count = REQ_count + 1
             totalsize = kmax * vmax
             rank      = Recv_info_p2r(I_prc_from,irank)
             tag       = Recv_info_p2r(I_prc_from,irank) + 2000000

             call MPI_ISEND( sendbuf_h2p_DP(1,I_SPL), & ! [IN]
                             totalsize,               & ! [IN]
                             MPI_DOUBLE_PRECISION,    & ! [IN]
                             rank,                    & ! [IN]
                             tag,                     & ! [IN]
                             PRC_LOCAL_COMM_WORLD,    & ! [IN]
                             REQ_list(REQ_count),     & ! [OUT]
                             ierr                     ) ! [OUT]
          endif
       enddo
    enddo

    !--- copy p2r-reverse
    do irank = 1, Copy_nmax_p2r
       do ipos = 1, Copy_info_p2r(I_size)
          ij_from = Copy_list_p2r(I_grid_to  ,ipos)
          l_from  = Copy_list_p2r(I_l_to     ,ipos)
          r_from  = PRC_RGN_lp2r(l_from,Copy_info_p2r(I_prc_to))
          ij_to   = Copy_list_p2r(I_grid_from,ipos)
          l_to    = Copy_list_p2r(I_l_from   ,ipos)
          r_to    = PRC_RGN_lp2r(l_to,Copy_info_p2r(I_prc_from))

          if (      r_from == PRC_RGN_rgn4pl(I_NPL) &
               .OR. r_from == PRC_RGN_rgn4pl(I_SPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                var_pl(ij_to,k,l_to,v) = var(ij_from,k,l_from,v)
             enddo
             enddo
          endif
       enddo
    enddo

    !--- wait all
    if ( REQ_count > 0 ) then
       call MPI_WAITALL( REQ_count,           & ! [IN]
                         REQ_list(:),         & ! [IN]
                         MPI_STATUSES_IGNORE, & ! [OUT]
                         ierr                 ) ! [OUT]
    endif

    !--- unpack p2r-reverse
    do irank = 1, Send_nmax_p2r
       do ipos = 1, Send_info_p2r(I_size,irank)
          l_from = Send_list_p2r(I_l_to     ,ipos,irank)
          r_from = PRC_RGN_lp2r(l_from,Send_info_p2r(I_prc_to,irank))

          ij_to  = Send_list_p2r(I_grid_from,ipos,irank)
          l_to   = Send_list_p2r(I_l_from   ,ipos,irank)

          if ( r_from == PRC_RGN_rgn4pl(I_NPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                kk = (v-1) * kmax + k
                var_pl(ij_to,k,l_to,v) = recvbuf_h2p_DP(kk,I_NPL)
             enddo
             enddo
          endif

          if ( r_from == PRC_RGN_rgn4pl(I_SPL) ) then
             do v = 1, vmax
             do k = 1, kmax
                kk = (v-1) * kmax + k
                var_pl(ij_to,k,l_to,v) = recvbuf_h2p_DP(kk,I_SPL)
             enddo
             enddo
          endif
       enddo
    enddo

    endif

    call COMM_data_transfer_DP(var,var_pl)

    call PROF_rapend('COMM_var',2)

    return
  end subroutine COMM_var_DP

  !-----------------------------------------------------------------------------
  !> suffix calculation
  !> @return suf
  function suf(i,j) result(suffix)
    implicit none

    integer  :: suffix
    integer  :: i, j
    !---------------------------------------------------------------------------

    suffix = ADM_gall_1d * (j-1) + i

  end function suf

  !-----------------------------------------------------------------------------
  subroutine COMM_debugtest
    use scale_prc, only: &
       PRC_myrank,   &
       PRC_MPIfinish
    use scale_prc_icoA, only: &
       PRC_RGN_l2r,      &
       PRC_RGN_vert_num, &
       PRC_RGN_have_sgp, &
       I_N,              &
       I_S
    implicit none

    real(RP) :: var   (ADM_gall   ,ADM_kall,ADM_lall   ,4)
    real(RP) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,4)

    integer  :: i, j, k, l, ij, rgnid, prc
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ TEST start'

    var   (:,:,:,:) = -999.0_RP
    var_pl(:,:,:,:) = -999.0_RP

    do l = 1, ADM_lall
       rgnid = PRC_RGN_l2r(l)
       prc   = PRC_myrank

       do k = ADM_kmin, ADM_kmax
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = ADM_gall_1d * (j-1) + i

             var(ij,k,l,1) = real(prc,  kind=RP)
             var(ij,k,l,2) = real(rgnid,kind=RP)
             var(ij,k,l,3) = real(i,    kind=RP)
             var(ij,k,l,4) = real(j,    kind=RP)
          enddo
          enddo
       enddo

       if ( PRC_RGN_have_sgp(l) ) then
          do k = ADM_kmin, ADM_kmax
             var(1,k,l,:) = -1.0_RP
          enddo
       endif
    enddo

    do l = 1, ADM_lall_pl
       rgnid = l
       prc   = PRC_myrank

       do k  = ADM_kmin, ADM_kmax
       do ij = 1, ADM_gall_pl
          var_pl(ij,k,l,1) = real(-prc,  kind=RP)
          var_pl(ij,k,l,2) = real(-rgnid,kind=RP)
          var_pl(ij,k,l,3) = real(-ij,   kind=RP)
          var_pl(ij,k,l,4) = real(-ij,   kind=RP)
       enddo
       enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '##### (prc,rgnid) #####'
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,1)),',',int(var(ij,k,l,2)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,1)),',',int(var_pl(ij,k,l,2)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '##### (i,j) #####'
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,3)),',',int(var(ij,k,l,4)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,3)),',',int(var_pl(ij,k,l,4)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo



    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Communication start'

    call COMM_data_transfer( var(:,:,:,:), var_pl(:,:,:,:) )

    if( IO_L ) write(IO_FID_LOG,*) '##### (prc,rgnid) #####'
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,1)),',',int(var(ij,k,l,2)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,1)),',',int(var_pl(ij,k,l,2)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '##### (i,j) #####'
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,3)),',',int(var(ij,k,l,4)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,3)),',',int(var_pl(ij,k,l,4)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo



    do l = 1, ADM_lall
       rgnid = PRC_RGN_l2r(l)
       prc   = PRC_myrank

       if ( PRC_RGN_vert_num(I_N,rgnid) == ADM_vlink ) then
          do k = ADM_kmin, ADM_kmax
             j  = ADM_gmax+1
             i  = ADM_gmin
             ij = ADM_gall_1d * (j-1) + i

             var(ij,k,l,1) = real(prc,  kind=RP)
             var(ij,k,l,2) = real(rgnid,kind=RP)
             var(ij,k,l,3) = real(i,    kind=RP)
             var(ij,k,l,4) = real(j,    kind=RP)
          enddo
       endif

       if ( PRC_RGN_vert_num(I_S,rgnid) == ADM_vlink ) then
          do k = ADM_kmin, ADM_kmax
             j  = ADM_gmin
             i  = ADM_gmax+1
             ij = ADM_gall_1d * (j-1) + i

             var(ij,k,l,1) = real(prc,  kind=RP)
             var(ij,k,l,2) = real(rgnid,kind=RP)
             var(ij,k,l,3) = real(i,    kind=RP)
             var(ij,k,l,4) = real(j,    kind=RP)
          enddo
       endif

    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ pole fill start'

    call COMM_var( var(:,:,:,:), var_pl(:,:,:,:), ADM_kall, 4 )

    if( IO_L ) write(IO_FID_LOG,*) '##### (prc,rgnid) #####'
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,1)),',',int(var(ij,k,l,2)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,1)),',',int(var_pl(ij,k,l,2)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '##### (i,j) #####'
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,3)),',',int(var(ij,k,l,4)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') '        |'
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, '|'
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,3)),',',int(var_pl(ij,k,l,4)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo

    call PRC_MPIfinish

    stop
  end subroutine COMM_debugtest

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_sum_SP( localsum, globalsum )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    real(SP), intent(in)  :: localsum
    real(SP), intent(out) :: globalsum

    real(SP) :: sendbuf(1)
    real(SP) :: recvbuf(PRC_nprocs)

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_pl ) then
       sendbuf(1) = localsum

       call MPI_Allgather( sendbuf,              &
                           1,                    &
                           MPI_REAL,             &
                           recvbuf,              &
                           1,                    &
                           MPI_REAL,             &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr                  )

       globalsum = sum( recvbuf(:) )
    else
       globalsum = localsum
    endif

    return
  end subroutine COMM_Stat_sum_SP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_sum_DP( localsum, globalsum )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    real(DP), intent(in)  :: localsum
    real(DP), intent(out) :: globalsum

    real(DP) :: sendbuf(1)
    real(DP) :: recvbuf(PRC_nprocs)

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_pl ) then
       sendbuf(1) = localsum

       call MPI_Allgather( sendbuf,              &
                           1,                    &
                           MPI_DOUBLE_PRECISION, &
                           recvbuf,              &
                           1,                    &
                           MPI_DOUBLE_PRECISION, &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr                  )

       globalsum = sum( recvbuf(:) )
    else
       globalsum = localsum
    endif

    return
  end subroutine COMM_Stat_sum_DP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_sum_eachlayer_SP( kall, localsum, globalsum )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    integer,  intent(in)  :: kall
    real(SP), intent(in)  :: localsum (kall)
    real(SP), intent(out) :: globalsum(kall)

    real(SP) :: sendbuf(kall)
    integer  :: displs (PRC_nprocs)
    integer  :: counts (PRC_nprocs)
    real(SP) :: recvbuf(kall,PRC_nprocs)

    integer  :: ierr
    integer  :: k, p
    !---------------------------------------------------------------------------

    do p = 1, PRC_nprocs
       displs(p) = (p-1) * kall
       counts(p) = kall
    enddo

    if ( COMM_pl ) then
       sendbuf(:) = localsum(:)

       call MPI_Allgatherv( sendbuf,              &
                            kall,                 &
                            MPI_REAL,             &
                            recvbuf,              &
                            counts,               &
                            displs,               &
                            MPI_REAL,             &
                            PRC_LOCAL_COMM_WORLD, &
                            ierr                  )

       do k = 1, kall
          globalsum(k) = sum( recvbuf(k,:) )
       enddo
    else
       do k = 1, kall
          globalsum(k) = localsum(k)
       enddo
    endif

    return
  end subroutine COMM_Stat_sum_eachlayer_SP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_sum_eachlayer_DP( kall, localsum, globalsum )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    integer,  intent(in)  :: kall
    real(DP), intent(in)  :: localsum (kall)
    real(DP), intent(out) :: globalsum(kall)

    real(DP) :: sendbuf(kall)
    integer  :: displs (PRC_nprocs)
    integer  :: counts (PRC_nprocs)
    real(DP) :: recvbuf(kall,PRC_nprocs)

    integer  :: ierr
    integer  :: k, p
    !---------------------------------------------------------------------------

    do p = 1, PRC_nprocs
       displs(p) = (p-1) * kall
       counts(p) = kall
    enddo

    if ( COMM_pl ) then
       sendbuf(:) = localsum(:)

       call MPI_Allgatherv( sendbuf,              &
                            kall,                 &
                            MPI_DOUBLE_PRECISION, &
                            recvbuf,              &
                            counts,               &
                            displs,               &
                            MPI_DOUBLE_PRECISION, &
                            PRC_LOCAL_COMM_WORLD, &
                            ierr                  )

       do k = 1, kall
          globalsum(k) = sum( recvbuf(k,:) )
       enddo
    else
       do k = 1, kall
          globalsum(k) = localsum(k)
       enddo
    endif

    return
  end subroutine COMM_Stat_sum_eachlayer_DP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_avg_SP( localavg, globalavg )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    real(SP), intent(in)  :: localavg
    real(SP), intent(out) :: globalavg

    real(SP) :: sendbuf(1)
    real(SP) :: recvbuf(PRC_nprocs)

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_pl ) then
       sendbuf(1) = localavg

       call MPI_Allgather( sendbuf,              &
                           1,                    &
                           MPI_REAL,             &
                           recvbuf,              &
                           1,                    &
                           MPI_REAL,             &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr                  )

       globalavg = sum( recvbuf(:) ) / real(PRC_nprocs,kind=SP)
    else
       globalavg = localavg
    endif

    return
  end subroutine COMM_Stat_avg_SP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_avg_DP( localavg, globalavg )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    real(DP), intent(in)  :: localavg
    real(DP), intent(out) :: globalavg

    real(DP) :: sendbuf(1)
    real(DP) :: recvbuf(PRC_nprocs)

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_pl ) then
       sendbuf(1) = localavg

       call MPI_Allgather( sendbuf,              &
                           1,                    &
                           MPI_DOUBLE_PRECISION, &
                           recvbuf,              &
                           1,                    &
                           MPI_DOUBLE_PRECISION, &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr                  )

       globalavg = sum( recvbuf(:) ) / real(PRC_nprocs,kind=DP)
    else
       globalavg = localavg
    endif

    return
  end subroutine COMM_Stat_avg_DP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_max_SP( localmax, globalmax )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    real(SP), intent(in)  :: localmax
    real(SP), intent(out) :: globalmax

    real(SP) :: sendbuf(1)
    real(SP) :: recvbuf(PRC_nprocs)

    integer  :: ierr
    !---------------------------------------------------------------------------

    sendbuf(1) = localmax

    call MPI_Allgather( sendbuf,              &
                        1,                    &
                        MPI_REAL,             &
                        recvbuf,              &
                        1,                    &
                        MPI_REAL,             &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )

    globalmax = maxval( recvbuf(:) )

    return
  end subroutine COMM_Stat_max_SP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_max_DP( localmax, globalmax )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    real(DP), intent(in)  :: localmax
    real(DP), intent(out) :: globalmax

    real(DP) :: sendbuf(1)
    real(DP) :: recvbuf(PRC_nprocs)

    integer  :: ierr
    !---------------------------------------------------------------------------

    sendbuf(1) = localmax

    call MPI_Allgather( sendbuf,              &
                        1,                    &
                        MPI_DOUBLE_PRECISION, &
                        recvbuf,              &
                        1,                    &
                        MPI_DOUBLE_PRECISION, &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )

    globalmax = maxval( recvbuf(:) )

    return
  end subroutine COMM_Stat_max_DP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_min_SP( localmin, globalmin )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    real(SP), intent(in)  :: localmin
    real(SP), intent(out) :: globalmin

    real(SP) :: sendbuf(1)
    real(SP) :: recvbuf(PRC_nprocs)

    integer  :: ierr
    !---------------------------------------------------------------------------

    sendbuf(1) = localmin

    call MPI_Allgather( sendbuf,        &
                        1,              &
                        MPI_REAL,       &
                        recvbuf,        &
                        1,              &
                        MPI_REAL,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr            )

    globalmin = minval( recvbuf(:) )

    return
  end subroutine COMM_Stat_min_SP

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_min_DP( localmin, globalmin )
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs
    implicit none

    real(DP), intent(in)  :: localmin
    real(DP), intent(out) :: globalmin

    real(DP) :: sendbuf(1)
    real(DP) :: recvbuf(PRC_nprocs)

    integer  :: ierr
    !---------------------------------------------------------------------------

    sendbuf(1) = localmin

    call MPI_Allgather( sendbuf,              &
                        1,                    &
                        MPI_DOUBLE_PRECISION, &
                        recvbuf,              &
                        1,                    &
                        MPI_DOUBLE_PRECISION, &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )

    globalmin = minval( recvbuf(:) )

    return
  end subroutine COMM_Stat_min_DP

end module scale_comm_icoA
