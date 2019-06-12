!-------------------------------------------------------------------------------
!> module process / icoA
!!
!! @par Description
!!          MPI process management module for Icosahedral-A grid
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_prc_icoA
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prc, only: &
     PRC_masterrank
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PRC_ICOA_setup
  public :: PRC_ICOA_RGN_generate

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !====== Information for processes ======

  integer,  public, parameter   :: PRC_rank_pl = PRC_masterrank !< process ID which manages the pole regions
  logical,  public              :: PRC_have_pl                  !< this ID manages pole region?

  !====== Information for processes-region relationship ======

  integer,  public, parameter   :: I_l   = 1                    !< local region
  integer,  public, parameter   :: I_prc = 2                    !< process

  integer,  public, parameter   :: I_RGNID = 1                  !< region id
  integer,  public, parameter   :: I_DIR   = 2                  !< direction

  ! Identifiers of directions of region edges
  integer,  public, parameter   :: I_SW = 1                     !< south west
  integer,  public, parameter   :: I_NW = 2                     !< north west
  integer,  public, parameter   :: I_NE = 3                     !< north east
  integer,  public, parameter   :: I_SE = 4                     !< south east

  ! Identifiers of directions of region vertices
  integer,  public, parameter   :: I_W = 1                      !< west
  integer,  public, parameter   :: I_N = 2                      !< north
  integer,  public, parameter   :: I_E = 3                      !< east
  integer,  public, parameter   :: I_S = 4                      !< south

  ! Identifier of poles (north pole or south pole)
  integer,  public, parameter   :: I_NPL = 1                    !< north pole
  integer,  public, parameter   :: I_SPL = 2                    !< south pole

  ! main parameter
  integer,  public              :: PRC_RGN_level    = -1        !< region division level
  integer,  public              :: PRC_RGN_ndiamond = 10        !< number of diamonds
  integer,  public              :: PRC_RGN_vlink    = 5         !< maximum number of vertex linkage, ICO:5

  ! region
  integer,  public              :: PRC_RGN_total                !< number of regular region (global total)
  integer,  public              :: PRC_RGN_local                !< number of regular region (local)
  integer,  public, parameter   :: PRC_RGN_total_pl  = 2        !< number of pole    region

  integer,  public, parameter   :: PRC_RGN_local_lim = 2560     !< maximum number of regular region (local)

  integer,  public, allocatable :: PRC_RGN_edge_tab   (:,:,:)   !< region link information (for 4 edges)

  integer,  public, allocatable :: PRC_RGN_vert_num   (:,:)     !< number of region around the vertex (4 vertexes)
  logical,  public, allocatable :: PRC_RGN_vert_pl    (:,:)     !< the northern/southern vertex is around the pole point?
  integer,  public, allocatable :: PRC_RGN_vert_tab   (:,:,:,:) !< region link information (for 4 vertexes)
  integer,  public, allocatable :: PRC_RGN_vert_tab_pl(:,:,:)   !< region link information (for 4 vertexes)

  integer,  public, allocatable :: PRC_RGN_lp2r(:,:)            !< l,prc       => rgn
  integer,  public, allocatable :: PRC_RGN_r2lp(:,:)            !< rgn         => l,prc
  integer,  public, allocatable :: PRC_RGN_l2r (:)              !< l,prc_me    => rgn

  integer,  public              :: PRC_RGN_r2p_pl(PRC_RGN_total_pl) !< process ID which have the pole regions
  integer,  public              :: PRC_RGN_rgn4pl(PRC_RGN_total_pl) !< region, having pole data in the halo

  !====== Information for regions ======

  logical,  public, allocatable :: PRC_RGN_have_sgp(:)          !< region have singlar point?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: PRC_ICOA_RGN_setup
  private :: PRC_ICOA_RGN_input
  private :: PRC_ICOA_RGN_output
  private :: PRC_ICOA_RGN_vertex_walkaround
  private :: output_info

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private            :: PRC_ICOA_RGN_in_fname  = '' !< input  file name for region management file
  character(len=H_LONG), private            :: PRC_ICOA_RGN_out_fname = '' !< output file name for region management file

  character(len=2),      private, parameter :: PRC_RGN_edgename(4) = (/'SW','NW','NE','SE'/) !< name table
  character(len=2),      private, parameter :: PRC_RGN_vertname(4) = (/'W ','N ','E ','S '/) !< name table

  logical, private :: debug = .false. !< debug option

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Processor topology
  subroutine PRC_ICOA_setup
    use scale_prc, only: &
       PRC_abort,                &
       PRC_ABORT_COMM_WORLD,     &
       PRC_UNIVERSAL_COMM_WORLD, &
       PRC_UNIVERSAL_nprocs,     &
       PRC_UNIVERSAL_myrank,     &
       PRC_UNIVERSAL_IsMaster,   &
       PRC_GLOBAL_COMM_WORLD,    &
       PRC_GLOBAL_nprocs,        &
       PRC_GLOBAL_myrank,        &
       PRC_GLOBAL_IsMaster,      &
       PRC_LOCAL_COMM_WORLD,     &
       PRC_nprocs,               &
       PRC_myrank,               &
       PRC_IsMaster
    implicit none

    namelist / PARAM_PRC_ICOA / &
       PRC_RGN_level,    &
       PRC_RGN_ndiamond, &
       debug

    integer :: l, rgnid

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("PRC_ICOA_setup",*) 'Setup'

    if ( IO_L ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start MPI'
       LOG_NEWLINE
       LOG_INFO("PRC_ICOA_setup",*) 'Process information '
       LOG_INFO_CONT('(1x,A,I12)')  'UNIVERSAL_COMM_WORLD        : ', PRC_UNIVERSAL_COMM_WORLD
       LOG_INFO_CONT('(1x,A,I12)')  'total process [UNIVERSAL]   : ', PRC_UNIVERSAL_nprocs
       LOG_INFO_CONT('(1x,A,I12)')  'my process ID [UNIVERSAL]   : ', PRC_UNIVERSAL_myrank
       LOG_INFO_CONT('(1x,A,L12)')  'master rank?  [UNIVERSAL]   : ', PRC_UNIVERSAL_IsMaster
       LOG_INFO_CONT('(1x,A,I12)')  'GLOBAL_COMM_WORLD           : ', PRC_GLOBAL_COMM_WORLD
       LOG_INFO_CONT('(1x,A,I12)')  'total process [GLOBAL]      : ', PRC_GLOBAL_nprocs
       LOG_INFO_CONT('(1x,A,I12)')  'my process ID [GLOBAL]      : ', PRC_GLOBAL_myrank
       LOG_INFO_CONT('(1x,A,L12)')  'master rank?  [GLOBAL]      : ', PRC_GLOBAL_IsMaster
       LOG_INFO_CONT('(1x,A,I12)')  'LOCAL_COMM_WORLD            : ', PRC_LOCAL_COMM_WORLD
       LOG_INFO_CONT('(1x,A,I12)')  'total process [LOCAL]       : ', PRC_nprocs
       LOG_INFO_CONT('(1x,A,I12)')  'my process ID [LOCAL]       : ', PRC_myrank
       LOG_INFO_CONT('(1x,A,L12)')  'master rank?  [LOCAL]       : ', PRC_IsMaster
       LOG_INFO_CONT('(1x,A,I12)')  'ABORT_COMM_WORLD            : ', PRC_ABORT_COMM_WORLD
       LOG_INFO_CONT('(1x,A,I12)')  'master rank ID [each world] : ', PRC_masterrank
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PRC_ICOA,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("PRC_ICOA_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("PRC_ICOA_setup",*) 'Not appropriate names in namelist PARAM_PRC_ICOA. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_PRC_ICOA)

    if ( PRC_RGN_level < 0 ) then
       LOG_ERROR("PRC_ICOA_setup",*) 'PRC_RGN_level is not appropriate :', PRC_RGN_level
       call PRC_abort
    endif

    if (      PRC_RGN_ndiamond ==  8 &
         .OR. PRC_RGN_ndiamond == 10 &
         .OR. PRC_RGN_ndiamond == 12 ) then
       PRC_RGN_vlink = PRC_RGN_ndiamond / 2
    else
       LOG_ERROR("PRC_ICOA_setup",*) 'PRC_RGN_ndiamond is not appropriate :', PRC_RGN_ndiamond
       call PRC_abort
    endif

    PRC_RGN_total = 2**PRC_RGN_level * 2**PRC_RGN_level * PRC_RGN_ndiamond
    PRC_RGN_local = PRC_RGN_total / PRC_nprocs

    if ( mod(PRC_RGN_total,PRC_nprocs) /= 0 ) then
       LOG_ERROR("PRC_ICOA_setup",*) 'Number of total region must be divisible by the number of process', PRC_RGN_total, PRC_nprocs
       call PRC_abort
    endif

    if ( PRC_RGN_local > PRC_RGN_local_lim ) then
       LOG_ERROR("PRC_ICOA_setup",*) 'limit exceed! local region: ', PRC_RGN_local, PRC_RGN_local_lim
       call PRC_abort
    endif

    call PRC_ICOA_RGN_setup

    ! pole region management flag
    if ( PRC_myrank == PRC_rank_pl ) then
       PRC_have_pl = .true.
    else
       PRC_have_pl = .false.
    endif

    ! singlar point management flag
    allocate( PRC_RGN_have_sgp(PRC_RGN_local) )
    PRC_RGN_have_sgp(:) = .false.

    do l = 1, PRC_RGN_local
       rgnid = PRC_RGN_lp2r(l,PRC_myrank)
       if ( PRC_RGN_vert_num(I_W,rgnid) == 3 ) then
          PRC_RGN_have_sgp(l) = .true.
       endif
    enddo

    call output_info

    return
  end subroutine PRC_ICOA_setup

  !-----------------------------------------------------------------------------
  subroutine PRC_ICOA_RGN_setup
    use scale_prc, only: &
       PRC_abort,  &
       PRC_nprocs, &
       PRC_myrank
    implicit none

    namelist / PARAM_PRC_ICOA_RGN / &
       PRC_ICOA_RGN_in_fname,  &
       PRC_ICOA_RGN_out_fname, &
       debug

    integer  :: l, p, r
    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("PRC_ICOA_RGN_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PRC_ICOA_RGN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("PRC_ICOA_RGN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("PRC_ICOA_RGN_setup",*) 'Not appropriate names in namelist PARAM_PRC_ICOA_RGN. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_PRC_ICOA_RGN)

    ! Global information (Each process has all the information)
    allocate( PRC_RGN_edge_tab(I_RGNID:I_DIR,I_SW:I_SE,PRC_RGN_total) )
    allocate( PRC_RGN_lp2r    (PRC_RGN_local,0:PRC_nprocs-1) )

    if ( PRC_ICOA_RGN_in_fname /= '' ) then
       call PRC_ICOA_RGN_input( PRC_ICOA_RGN_in_fname,   & ! [IN]
                                PRC_nprocs,              & ! [IN]
                                PRC_RGN_total,           & ! [IN]
                                PRC_RGN_local,           & ! [IN]
                                PRC_RGN_edge_tab(:,:,:), & ! [OUT]
                                PRC_RGN_lp2r    (:,:)    ) ! [OUT]
    else
       LOG_INFO("PRC_ICOA_RGN_setup",*) 'input file is not specified.'

       call PRC_ICOA_RGN_generate( PRC_RGN_level,           & ! [IN]
                                   PRC_RGN_ndiamond,        & ! [IN]
                                   PRC_nprocs,              & ! [IN]
                                   PRC_RGN_total,           & ! [IN]
                                   PRC_RGN_local,           & ! [IN]
                                   PRC_RGN_edge_tab(:,:,:), & ! [OUT]
                                   PRC_RGN_lp2r    (:,:)    ) ! [OUT]
    endif

    if ( PRC_ICOA_RGN_out_fname /= '' ) then
       call PRC_ICOA_RGN_output( PRC_ICOA_RGN_out_fname,  & ! [IN]
                                 PRC_nprocs,              & ! [IN]
                                 PRC_RGN_total,           & ! [IN]
                                 PRC_RGN_local,           & ! [IN]
                                 PRC_RGN_edge_tab(:,:,:), & ! [IN]
                                 PRC_RGN_lp2r    (:,:)    ) ! [IN]
    endif

    !--- additional table (reversal,local)
    allocate( PRC_RGN_r2lp(2,PRC_RGN_total) )
    allocate( PRC_RGN_l2r (  PRC_RGN_local) )

    do p = 0, PRC_nprocs-1
    do l = 1, PRC_RGN_local
       PRC_RGN_r2lp(I_l,  PRC_RGN_lp2r(l,p)) = l
       PRC_RGN_r2lp(I_prc,PRC_RGN_lp2r(l,p)) = p
    enddo
    enddo

    do l = 1, PRC_RGN_local
       PRC_RGN_l2r(l) = PRC_RGN_lp2r(l,PRC_myrank)
    enddo

    !--- region connection chains around the diamond vertexes
    allocate( PRC_RGN_vert_num   (I_W:I_S,PRC_RGN_total) )
    allocate( PRC_RGN_vert_pl    (PRC_RGN_total_pl,PRC_RGN_total) )
    allocate( PRC_RGN_vert_tab   (I_RGNID:I_DIR,I_W:I_S,PRC_RGN_total   ,PRC_RGN_vlink) )
    allocate( PRC_RGN_vert_tab_pl(I_RGNID:I_DIR,        PRC_RGN_total_pl,PRC_RGN_vlink) )

    call PRC_ICOA_RGN_vertex_walkaround( PRC_RGN_total,             & ! [IN]
                                         PRC_RGN_total_pl,          & ! [IN]
                                         PRC_RGN_vlink,             & ! [IN]
                                         PRC_RGN_edge_tab(:,:,:),   & ! [IN]
                                         PRC_RGN_vert_num(:,:),     & ! [OUT]
                                         PRC_RGN_vert_pl(:,:),      & ! [OUT]
                                         PRC_RGN_vert_tab(:,:,:,:), & ! [OUT]
                                         PRC_RGN_vert_tab_pl(:,:,:) ) ! [OUT]

    !--- tables for pole
    do r = 1, PRC_RGN_total
       if ( PRC_RGN_vert_pl(I_NPL,r) ) then
          PRC_RGN_rgn4pl(I_NPL) = r
          exit
       endif
    enddo

    do r = 1, PRC_RGN_total
       if ( PRC_RGN_vert_pl(I_SPL,r) ) then
          PRC_RGN_rgn4pl(I_SPL) = r
          exit
       endif
    enddo

    PRC_RGN_r2p_pl(I_NPL) = PRC_rank_pl
    PRC_RGN_r2p_pl(I_SPL) = PRC_rank_pl

    return
  end subroutine PRC_ICOA_RGN_setup

  !-----------------------------------------------------------------------------
  !> Input mnginfo file
  subroutine PRC_ICOA_RGN_input( &
       in_fname, &
       pall,     &
       rall,     &
       lall,     &
       edge_tab, &
       lp2r      )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in)  :: in_fname            !< input file
    integer,          intent(in)  :: pall                !< number of process        (global total)
    integer,          intent(in)  :: rall                !< number of regular region (global total)
    integer,          intent(in)  :: lall                !< number of regular region (local)
    integer,          intent(out) :: edge_tab(2,4,rall)  !< region link information (for 4 edges)
    integer,          intent(out) :: lp2r(lall,0:pall-1) !< l,prc => region

    integer :: num_of_rgn !< number of region

    namelist / rgn_info / &
       num_of_rgn

    integer :: rgnid      !< region ID
    integer :: sw(2) = -1 !< south-west region info
    integer :: nw(2) = -1 !< north-west region info
    integer :: ne(2) = -1 !< north-east region info
    integer :: se(2) = -1 !< south-east region info

    namelist / rgn_link_info / &
       rgnid, &
       sw,    &
       nw,    &
       ne,    &
       se

    integer :: num_of_proc !< number of processes

    namelist / proc_info / &
       num_of_proc

    integer :: peid                         !< process ID
    integer :: num_of_mng                   !< number of regions be managed
    integer :: mng_rgnid(PRC_RGN_local_lim) !< managed region ID

    namelist / rgn_mng_info / &
       peid,       &
       num_of_mng, &
       mng_rgnid

    integer  :: fid, ierr
    integer  :: r, p, l
    !---------------------------------------------------------------------------

    LOG_INFO("PRC_ICOA_RGN_input",*) 'input region management information file: ', trim(in_fname)

    fid = IO_get_available_fid()
    open( unit   = fid,            &
          file   = trim(in_fname), &
          form   = 'formatted',    &
          status = 'old',          &
          iostat = ierr            )

       ! ERROR if filename are not defined
       if ( ierr /= 0 ) then
          LOG_ERROR("PRC_ICOA_RGN_input",*) 'File is not found!', trim(in_fname)
          call PRC_abort
       endif

       read(fid,nml=rgn_info)

       if ( num_of_rgn /= rall ) then
          LOG_ERROR("PRC_ICOA_RGN_input",*) 'Missmatch of region number!'
          LOG_ERROR_CONT(*)                 'rall= ', rall,', num_of_rgn=', num_of_rgn
          call PRC_abort
       endif

       do r = 1, rall
          read(fid,nml=rgn_link_info)

          edge_tab(I_RGNID:I_DIR,I_SW,rgnid) = sw(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_NW,rgnid) = nw(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_NE,rgnid) = ne(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_SE,rgnid) = se(I_RGNID:I_DIR)
       enddo

       read(fid,nml=proc_info)
       if ( num_of_proc /= pall ) then
          LOG_ERROR("PRC_ICOA_RGN_input",*) ' Missmatch of process number!'
          LOG_ERROR_CONT(*)                 ' pall= ', pall, ', num_of_proc=', num_of_proc
          call PRC_abort
       endif

       do p = 0, pall-1
          read(fid,nml=rgn_mng_info)

          if ( p /= peid ) then
             LOG_ERROR("PRC_ICOA_RGN_input",*) 'Wrong peid: ', p, peid
             call PRC_abort
          endif

          if ( num_of_mng /= lall ) then
             LOG_ERROR("PRC_ICOA_RGN_input",*) 'number of local region is not match: ', p, num_of_mng, lall
             call PRC_abort
          endif

          lp2r(:,p) = -1 ! initialize
          do l = 1, lall
             lp2r(l,p) = mng_rgnid(l)
          enddo
       enddo

    close(fid)

    return
  end subroutine PRC_ICOA_RGN_input

  !-----------------------------------------------------------------------------
  !> Output mnginfo file
  subroutine PRC_ICOA_RGN_output( &
       out_fname, &
       pall,      &
       rall,      &
       lall,      &
       edge_tab,  &
       lp2r       )
    use scale_prc, only: &
       PRC_IsMaster
    implicit none

    character(len=*), intent(in) :: out_fname          !< output file
    integer,          intent(in) :: pall               !< number of process        (global total)
    integer,          intent(in) :: rall               !< number of regular region (global total)
    integer,          intent(in) :: lall               !< number of regular region (local)
    integer,          intent(in) :: edge_tab(2,4,rall) !< region link information (for 4 edges)
    integer,          intent(in) :: lp2r(lall,pall)    !< l,prc => region

    integer :: num_of_rgn !< number of region

    namelist / rgn_info / &
       num_of_rgn

    integer :: rgnid      !< region ID
    integer :: sw(2) = -1 !< south-west region info
    integer :: nw(2) = -1 !< north-west region info
    integer :: ne(2) = -1 !< north-east region info
    integer :: se(2) = -1 !< south-east region info

    namelist / rgn_link_info / &
       rgnid, &
       sw,    &
       nw,    &
       ne,    &
       se

    integer :: num_of_proc !< number of processes

    namelist / proc_info / &
       num_of_proc

    integer :: peid                         !< process ID
    integer :: num_of_mng                   !< number of regions be managed
    integer :: mng_rgnid(PRC_RGN_local_lim) !< managed region ID

    namelist / rgn_mng_info / &
       peid,       &
       num_of_mng, &
       mng_rgnid

    integer  :: fid
    integer  :: r, p, l
    !---------------------------------------------------------------------------

    if ( PRC_IsMaster ) then
       LOG_INFO("PRC_ICOA_RGN_output",*) 'output region management information file: ', trim(out_fname)

       fid = IO_get_available_fid()
       open( unit   = fid,             &
             file   = trim(out_fname), &
             form   = 'formatted'      )

          num_of_rgn = rall
          write(fid,nml=rgn_info)

          do r = 1, rall
             rgnid = r
             sw(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_SW,rgnid)
             nw(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_NW,rgnid)
             ne(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_NE,rgnid)
             se(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_SE,rgnid)

             write(fid,nml=rgn_link_info)
          enddo

          num_of_proc = pall
          write(fid,nml=proc_info)

          do p = 0, pall-1
             peid       = p
             num_of_mng = lall

             mng_rgnid(:) = -1
             do l = 1, lall
                mng_rgnid(l) = lp2r(l,p)
             enddo

             write(fid,nml=rgn_mng_info)
          enddo

       close(fid)
    else
       LOG_INFO("PRC_ICOA_RGN_output",*) 'output region management information file at the master process'
    endif

    return
  end subroutine PRC_ICOA_RGN_output

  !-----------------------------------------------------------------------------
  !> Generate region management info
  Subroutine PRC_ICOA_RGN_generate( &
       rlevel,   &
       ndmd,     &
       pall,     &
       rall,     &
       lall,     &
       edge_tab, &
       lp2r      )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in)  :: rlevel              !< region division level
    integer, intent(in)  :: ndmd                !< number of diamonds
    integer, intent(in)  :: pall                !< number of process        (global total)
    integer, intent(in)  :: rall                !< number of regular region (global total)
    integer, intent(in)  :: lall                !< number of regular region (local)
    integer, intent(out) :: edge_tab(2,4,rall)  !< region link information (for 4 edges)
    integer, intent(out) :: lp2r(lall,0:pall-1) !< l,prc => region

    integer :: dmd_data(4,ndmd)
    integer :: rall_1d, rall_1dmd

    integer :: d_nb, i_nb, j_nb, rgnid_nb, direction
    integer :: d, i, j, rgnid
    integer :: l, p
    !---------------------------------------------------------------------------

    LOG_INFO("PRC_ICOA_RGN_generate",*) 'generate region management information file'

    if    ( ndmd == 10 ) then
       LOG_INFO_CONT(*) 'Topology: icosahedral'
       dmd_data(:, 1) = (/  6, 5, 2,10 /)
       dmd_data(:, 2) = (/ 10, 1, 3, 9 /)
       dmd_data(:, 3) = (/  9, 2, 4, 8 /)
       dmd_data(:, 4) = (/  8, 3, 5, 7 /)
       dmd_data(:, 5) = (/  7, 4, 1, 6 /)
       dmd_data(:, 6) = (/  7, 5, 1,10 /)
       dmd_data(:, 7) = (/  8, 4, 5, 6 /)
       dmd_data(:, 8) = (/  9, 3, 4, 7 /)
       dmd_data(:, 9) = (/ 10, 2, 3, 8 /)
       dmd_data(:,10) = (/  6, 1, 2, 9 /)
    elseif( ndmd == 12 ) then
       LOG_INFO_CONT(*) 'Topology: icosatetrahedral'
       dmd_data(:, 1) = (/  7, 6, 2,12 /)
       dmd_data(:, 2) = (/ 12, 1, 3,11 /)
       dmd_data(:, 3) = (/ 11, 2, 4,10 /)
       dmd_data(:, 4) = (/ 10, 3, 5, 9 /)
       dmd_data(:, 5) = (/  9, 4, 6, 8 /)
       dmd_data(:, 6) = (/  8, 5, 1, 7 /)
       dmd_data(:, 7) = (/  8, 6, 1,12 /)
       dmd_data(:, 8) = (/  9, 5, 6, 7 /)
       dmd_data(:, 9) = (/ 10, 4, 5, 8 /)
       dmd_data(:,10) = (/ 11, 3, 4, 9 /)
       dmd_data(:,11) = (/ 12, 2, 3,10 /)
       dmd_data(:,12) = (/  7, 1, 2,11 /)
    endif

    rall_1d   = 2**rlevel
    rall_1dmd = rall_1d*rall_1d

    !--- make region link table
    do d = 1, ndmd
    do i = 1, rall_1d
    do j = 1, rall_1d
       rgnid = (d-1)*rall_1dmd + (j-1)*rall_1d + i

       !--- I_SW
       if ( j == 1 ) then
          if ( d <= ndmd / 2 ) then
             i_nb = i
             j_nb = rall_1d
             d_nb = dmd_data(I_SW,d)
             direction = I_NE
          else
             i_nb = rall_1d
             j_nb = rall_1d+1-i
             d_nb = dmd_data(I_SW,d)
             direction = I_SE
          endif
       else
          i_nb = i
          j_nb = j-1
          d_nb = d
          direction = I_NE
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_SW,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_SW,rgnid) = direction

       !--- I_NW
       if ( i == 1 ) then
          if ( d <= ndmd / 2 ) then
             i_nb = rall_1d+1-j
             j_nb = rall_1d
             d_nb = dmd_data(I_NW,d)
             direction = I_NE
          else
             i_nb = rall_1d
             j_nb = j
             d_nb = dmd_data(I_NW,d)
             direction = I_SE
          endif
       else
          i_nb = i-1
          j_nb = j
          d_nb = d
          direction = I_SE
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_NW,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_NW,rgnid) = direction

       !--- I_NE
       if ( j == rall_1d ) then
          if ( d <= ndmd / 2 ) then
             i_nb = 1
             j_nb = rall_1d+1-i
             d_nb = dmd_data(I_NE,d)
             direction = I_NW
          else
             i_nb = i
             j_nb = 1
             d_nb = dmd_data(I_NE,d)
             direction = I_SW
          endif
       else
          i_nb = i
          j_nb = j+1
          d_nb = d
          direction = I_SW
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_NE,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_NE,rgnid) = direction

       !--- I_SE
       if ( i == rall_1d ) then
          if ( d <= ndmd / 2 ) then
             i_nb = 1
             j_nb = j
             d_nb = dmd_data(I_SE,d)
             direction = I_NW
          else
             i_nb = rall_1d+1-j
             j_nb = 1
             d_nb = dmd_data(I_SE,d)
             direction = I_SW
          endif
       else
          i_nb = i+1
          j_nb = j
          d_nb = d
          direction = I_NW
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_SE,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_SE,rgnid) = direction
    enddo
    enddo
    enddo

    !--- make region-pe relationship
    lp2r(:,:) = -1

    rgnid = 0
    do p = 0, pall-1
    do l = 1, lall
       rgnid = rgnid + 1

       lp2r(l,p) = rgnid
    enddo
    enddo

    return
  end Subroutine PRC_ICOA_RGN_generate

  !-----------------------------------------------------------------------------
  !> Search region with walking around the diamond vertexes
  subroutine PRC_ICOA_RGN_vertex_walkaround( &
       rall,       &
       rall_pl,    &
       vlink,      &
       edge_tab,   &
       vert_num,   &
       vert_pl,    &
       vert_tab,   &
       vert_tab_pl )
    implicit none

    integer, intent(in)  :: rall
    integer, intent(in)  :: rall_pl
    integer, intent(in)  :: vlink
    integer, intent(in)  :: edge_tab(2,4,rall)

    integer, intent(out) :: vert_num   (4,rall)
    logical, intent(out) :: vert_pl    (rall_pl,rall)
    integer, intent(out) :: vert_tab   (2,4,rall   ,vlink)
    integer, intent(out) :: vert_tab_pl(2  ,rall_pl,vlink)

    integer  :: rgnid, dir
    integer  :: rgnid_next, dir_next
    logical  :: IsAroundPole

    integer  :: r, d, v
    !---------------------------------------------------------------------------

    vert_num(:,:)     = -1
    vert_tab(:,:,:,:) = -1

    do r = 1, rall
    do d = I_W, I_S

       rgnid = r
       select case(d)
       case(I_W)
          dir = I_SW
       case(I_N)
          dir = I_NW
       case(I_E)
          dir = I_NE
       case(I_S)
          dir = I_SE
       end select

       v = 0
       do
          rgnid_next = edge_tab(I_RGNID,dir,rgnid)
          dir_next   = edge_tab(I_DIR,  dir,rgnid) - 1

          if( dir_next == 0 ) dir_next = 4
          v = v + 1
          vert_tab(I_RGNID,d,r,v) = rgnid
          vert_tab(I_DIR,  d,r,v) = dir

          rgnid = rgnid_next
          dir   = dir_next

          if( rgnid == r ) exit
       enddo
       vert_num(d,r) = v

    enddo
    enddo

    vert_pl    (:,:)   = .false.

    do r = 1, rall
       if ( vert_num(I_N,r) == vlink ) then
          IsAroundPole = .true.
          do v = 1, vlink
             if( vert_tab(I_DIR,I_N,r,v) /= I_N ) IsAroundPole = .false.
          enddo

          if ( IsAroundPole ) then
             vert_pl(I_NPL,r) = .true.
          endif
       endif

       if ( vert_num(I_S,r) == vlink ) then
          IsAroundPole = .true.
          do v = 1, vlink
             if( vert_tab(I_DIR,I_S,r,v) /= I_S ) IsAroundPole = .false.
          enddo

          if ( IsAroundPole ) then
             vert_pl(I_SPL,r) = .true.
          endif
       endif
    enddo

    vert_tab_pl(:,:,:) = -1

    do r = 1, rall
       if ( vert_pl(I_NPL,r) ) then
          do v = 1, vlink
             vert_tab_pl(I_RGNID,I_NPL,v) = vert_tab(I_RGNID,I_N,r,v)
             vert_tab_pl(I_DIR,  I_NPL,v) = vert_tab(I_DIR,  I_N,r,v)
          enddo
          exit
       endif
    enddo

    do r = 1, rall
       if ( vert_pl(I_SPL,r) ) then
          do v = 1, vlink
             vert_tab_pl(I_RGNID,I_SPL,v) = vert_tab(I_RGNID,I_S,r,v)
             vert_tab_pl(I_DIR,  I_SPL,v) = vert_tab(I_DIR,  I_S,r,v)
          enddo
          exit
       endif
    enddo

    return
  end subroutine PRC_ICOA_RGN_vertex_walkaround

  !-----------------------------------------------------------------------------
  subroutine output_info
    use scale_prc, only: &
       PRC_myrank
    implicit none

    integer          :: rgnid, rgnid_next
    character(len=2) :: dstr, dstr_next

    integer  :: l, d, v
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("PRC_ICOA_RGN_setup",'(1x,A)') 'Region management information'
    LOG_INFO_CONT('(1x,A,A)' )              'Grid sysytem                      : Icosahedral'
    LOG_INFO_CONT('(1x,A,I7)')              'number of diamond                 : ', PRC_RGN_ndiamond
    LOG_INFO_CONT('(1x,A,I7)')              'maximum number of vertex linkage  : ', PRC_RGN_vlink
    LOG_NEWLINE
    LOG_INFO_CONT('(1x,A,I7)')              'Region division level (RL)        : ', PRC_RGN_level
    LOG_INFO_CONT('(1x,A,I7,3(A,I4),A)')    'Total number of regular region    : ', PRC_RGN_total, &
                                            ' (', 2**PRC_RGN_level, ' x', 2**PRC_RGN_level, ' x', PRC_RGN_ndiamond, ' )'
    LOG_INFO_CONT('(1x,A,I7)')              '#  of region per process          : ', PRC_RGN_local
    LOG_INFO_CONT('(1x,A)'   )              'ID of region in my process        : '
    LOG_INFO_CONT(*)                        PRC_RGN_lp2r(:,PRC_myrank)
    LOG_INFO_CONT('(1x,A,I7)')              'Region ID, containing north pole  : ', PRC_RGN_rgn4pl(I_NPL)
    LOG_INFO_CONT('(1x,A,I7)')              'Region ID, containing south pole  : ', PRC_RGN_rgn4pl(I_SPL)
    LOG_INFO_CONT('(1x,A,I7)')              'Process rank, managing north pole : ', PRC_RGN_r2p_pl(I_NPL)
    LOG_INFO_CONT('(1x,A,I7)')              'Process rank, managing south pole : ', PRC_RGN_r2p_pl(I_SPL)

    if ( debug ) then
       LOG_NEWLINE
       LOG_INFO("PRC_ICOA_RGN_setup",'(1x,A)') 'Detailed region management information'

       LOG_INFO_CONT(*)'--- (l,myrank) => (rgn)'
       do l = 1, PRC_RGN_local
          rgnid = PRC_RGN_l2r(l)
          LOG_INFO_CONT('(1x,A,I4,A,I6,A,I6,A)') '--- (',l,',',PRC_myrank,') => (',rgnid,') '
       enddo

       LOG_NEWLINE
       LOG_INFO_CONT(*)'--- Link information'
       do l = 1, PRC_RGN_local
          rgnid = PRC_RGN_l2r(l)

          LOG_NEWLINE
          LOG_INFO_CONT(*)'--- edge link: (rgn,direction)'
          do d = I_SW, I_SE
             rgnid_next = PRC_RGN_edge_tab(I_RGNID,d,rgnid)
             dstr       = PRC_RGN_edgename(d)
             dstr_next  = PRC_RGN_edgename(PRC_RGN_edge_tab(I_DIR,d,rgnid))
             LOG_INFO_CONT('(5x,A,I6,A,A,A,I6,A,A,A)') '(',rgnid,',',dstr,') -> (', rgnid_next,',', dstr_next,')'
          enddo

          LOG_INFO_CONT(*)'--- vertex link: (rgn)'
          do d = I_W, I_S
             dstr = PRC_RGN_vertname(d)
             LOG_INFO_CONTNA('(5x,A,I6,A,A,A)') '(',rgnid,',',dstr,')'
             do v = 2, PRC_RGN_vert_num(d,rgnid)
                dstr = PRC_RGN_vertname(PRC_RGN_vert_tab(I_DIR,d,rgnid,v))
                LOG_INFO_CONTNA('(A,I6,A,A,A)') ' -> (',PRC_RGN_vert_tab(I_RGNID,d,rgnid,v),',',dstr,')'
             enddo
             LOG_NEWLINE
          enddo
       enddo

       LOG_NEWLINE
       LOG_INFO_CONT(*)'--- Pole information (in the global scope)'

       LOG_NEWLINE
       LOG_INFO_CONT(*)'--- Region ID, containing north pole data : ', PRC_RGN_rgn4pl(I_NPL)
       LOG_INFO_CONT(*)'--- vertex link: (north pole)'
       do v = 2, PRC_RGN_vlink
          rgnid = PRC_RGN_vert_tab_pl(I_RGNID,I_NPL,v)
          dstr  = PRC_RGN_vertname(PRC_RGN_vert_tab_pl(I_DIR,I_NPL,v))
          LOG_INFO_CONTNA('(A,I6,A,A,A)') ' -> (',rgnid,',',dstr,')'
       enddo
       rgnid = PRC_RGN_vert_tab_pl(I_RGNID,I_NPL,1)
       dstr  = PRC_RGN_vertname(PRC_RGN_vert_tab_pl(I_DIR,I_NPL,1))
       LOG_INFO_CONTNA('(A,I6,A,A,A)') ' -> (',rgnid,',',dstr,')'
       LOG_NEWLINE
       LOG_INFO_CONT(*)'--- process, managing north pole : ', PRC_RGN_r2p_pl(I_NPL)

       LOG_NEWLINE
       LOG_INFO_CONT(*)'--- Region ID, containing south pole data : ', PRC_RGN_rgn4pl(I_SPL)
       LOG_INFO_CONT(*)'--- vertex link: (south pole)'
       do v = 2, PRC_RGN_vlink
          rgnid = PRC_RGN_vert_tab_pl(I_RGNID,I_SPL,v)
          dstr  = PRC_RGN_vertname(PRC_RGN_vert_tab_pl(I_DIR,I_SPL,v))
          LOG_INFO_CONTNA('(A,I6,A,A,A)') ' -> (',rgnid,',',dstr,')'
       enddo
       rgnid = PRC_RGN_vert_tab_pl(I_RGNID,I_SPL,1)
       dstr  = PRC_RGN_vertname(PRC_RGN_vert_tab_pl(I_DIR,I_SPL,1))
       LOG_INFO_CONTNA('(A,I6,A,A,A)') ' -> (',rgnid,',',dstr,')'
       LOG_NEWLINE
       LOG_INFO_CONT(*)'--- process, managing south pole : ', PRC_RGN_r2p_pl(I_SPL)
    endif

    return
  end subroutine output_info

end module scale_prc_icoA
