!-------------------------------------------------------------------------------
!> Module administration
!!
!! @par Description
!!         This module is for the management of process and region
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_adm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ADM_proc_init
  public :: ADM_proc_stop
  public :: ADM_proc_finish
  public :: ADM_setup
  public :: ADM_mk_suffix
  public :: ADM_MPItime

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  ! Character length of system control
  integer, public, parameter :: ADM_NSYS = 32
  !
  ! Maximum length of file name
  integer, public, parameter :: ADM_MAXFNAME = 128

  !
  !====== Basic definition & information ======
  !
  ! Log file ID & Control file ID
  integer, public      :: ADM_LOG_FID = 6 ! default is STDOUT
  integer, public      :: ADM_CTL_FID = 35
  !
  ! Identifier for single computation or parallel computation
  integer, public, parameter :: ADM_MULTI_PRC  = 1
  !
  ! Identifiers of directions of region edges
  integer, public, parameter :: ADM_SW = 1
  integer, public, parameter :: ADM_NW = 2
  integer, public, parameter :: ADM_NE = 3
  integer, public, parameter :: ADM_SE = 4
  !
  ! Identifiers of directions of region vertices
  integer, public, parameter :: ADM_W = 1
  integer, public, parameter :: ADM_N = 2
  integer, public, parameter :: ADM_E = 3
  integer, public, parameter :: ADM_S = 4
  !
  !--- Identifier of triangle element (i-axis-side or j-axis side)
  integer, public, parameter :: ADM_TI = 1
  integer, public, parameter :: ADM_TJ = 2
  !
  !--- Identifier of line element (i-axis-side, ij-axis side, or j-axis side)
  integer, public, parameter :: ADM_AI  = 1
  integer, public, parameter :: ADM_AIJ = 2
  integer, public, parameter :: ADM_AJ  = 3
  !
  ! Identifier of 1 variable
  integer, public, parameter :: ADM_KNONE = 1

  ! Identifier of poles (north pole or south pole)
  integer, public, parameter :: ADM_NPL = 1
  integer, public, parameter :: ADM_SPL = 2

  ! Fist colomn on the table for region and direction
  integer, public, parameter :: ADM_RID = 1
  integer, public, parameter :: ADM_DIR = 2

#ifdef _FIXEDINDEX_
  include "inc_index.h"
#else
  !#############################################################################
  ! Basic Index Parameters
  !#############################################################################

  ! main parameter
  integer, public            :: ADM_glevel          ! grid   division level
  integer, public            :: ADM_rlevel          ! region division level
  integer, public            :: ADM_vlayer          ! number of vertical layer
  integer, public            :: ADM_DMD             ! number of diamond
  integer, public            :: ADM_prc_all         ! number of MPI process

  ! region
  integer, public            :: ADM_rgn_nmax        ! number of regular region
  integer, public            :: ADM_lall            ! number of regular region per process
  integer, public, parameter :: ADM_rgn_nmax_pl = 2 ! number of pole    region
  integer, public, parameter :: ADM_lall_pl     = 2 ! number of pole    region per process

  ! horizontal grid
  integer, public            :: ADM_gall            ! number of horizontal grid per regular region
  integer, public            :: ADM_gall_in         ! number of horizontal grid (inner part)
  integer, public            :: ADM_gall_1d         ! number of horizontal grid (1D)
  integer, public            :: ADM_gmin            ! start index of 1D horizontal grid
  integer, public            :: ADM_gmax            ! end   index of 1D horizontal grid

  integer, public            :: ADM_gall_pl         ! number of horizontal grid for pole region
  integer, public, parameter :: ADM_gslf_pl     = 1 ! index for pole point
  integer, public, parameter :: ADM_gmin_pl     = 2 ! start index of grid around the pole point
  integer, public            :: ADM_gmax_pl         ! end   index of grid around the pole point

  ! vertical grid
  integer, public            :: ADM_kall            ! number of vertical grid
  integer, public            :: ADM_kmin            ! start index of vertical grid
  integer, public            :: ADM_kmax            ! end   index of vertical grid

  ! List vectors
  integer, public            :: ADM_IooJoo_nmax
  integer, public            :: ADM_IooJmo_nmax
  integer, public            :: ADM_IooJop_nmax
  integer, public            :: ADM_IooJmp_nmax
  integer, public            :: ADM_ImoJoo_nmax
  integer, public            :: ADM_ImoJmo_nmax
  integer, public            :: ADM_ImoJop_nmax
  integer, public            :: ADM_ImoJmp_nmax
  integer, public            :: ADM_IopJoo_nmax
  integer, public            :: ADM_IopJmo_nmax
  integer, public            :: ADM_IopJop_nmax
  integer, public            :: ADM_IopJmp_nmax
  integer, public            :: ADM_ImpJoo_nmax
  integer, public            :: ADM_ImpJmo_nmax
  integer, public            :: ADM_ImpJop_nmax
  integer, public            :: ADM_ImpJmp_nmax
#endif

  !
  !====== Information for processes ======
  !
  integer, public            :: ADM_COMM_WORLD          ! communication world per member
  logical, public            :: ADM_MPI_alive = .false. ! MPI is alive?

  integer, public            :: ADM_prc_me              ! my process ID
  integer, public, parameter :: ADM_prc_run_master = 1  ! master process ID
  integer, public            :: ADM_prc_pl              ! process ID which manages the pole regions
  logical, public            :: ADM_have_pl             ! this ID manages pole region?

  integer, public            :: ADM_prc_npl             ! process ID which have the pole regions
  integer, public            :: ADM_prc_spl             ! process ID which have the pole regions
  integer, public            :: ADM_prc_nspl(ADM_NPL:ADM_SPL)

  !
  !====== Information for processes-region relationship ======
  !
  character(len=ADM_MAXFNAME), public :: ADM_rgnmngfname ! file name for region management info

  integer, public, parameter   :: PRC_RGN_NMAX   = 2560  ! maximum number of region per process.
  integer, public              :: ADM_vlink_nmax = -1    ! maximum number of vertex linkage
                                                         ! [XTMS] ICO:5, PSP:6, LCP, MLCP:k

  integer, public, allocatable :: ADM_prc_rnum(:)        ! number of regions managed by each process = ADM_lall
  integer, public, allocatable :: ADM_prc_tab (:,:)      ! table  of regions managed by each process

  integer, public, allocatable :: ADM_rgn2prc (:)        ! Table of process ID from region ID
                                                         ! ADM_rgn2prc(ADM_rgn_nmax)
  integer, public, allocatable :: ADM_rgn_etab(:,:,:)    ! table  of edge link information
                                                         ! ADM_rgn_etab( ADM_RID:ADM_DIR,ADM_SW:ADM_SE,ADM_rgn_nmax )
  integer, public, allocatable :: ADM_rgn_vnum(:,:)      ! Table of n-vertex-link at the region vertex
                                                         ! ADM_rgn_vnum( ADM_W:ADM_S,ADM_rgn_nmax )
  integer, public, allocatable :: ADM_rgn_vtab(:,:,:,:)  ! Table of vertex link information
                                                         ! ADM_rgn_vtab   ( ADM_RID:ADM_DIR, &
                                                         !                  ADM_W:ADM_S,     &
                                                         !                  ADM_rgn_nmax,    &
                                                         !                  ADM_vlink_nmax   )
  integer, public, allocatable :: ADM_rgn_vtab_pl(:,:,:) ! Table of vertex link information for poles
                                                         ! ADM_rgn_vtab_pl( ADM_RID:ADM_DIR, &
                                                         !                  ADM_rgn_nmax_pl, &
                                                         !                  ADM_vlink_nmax   )

  integer, public              :: ADM_rgnid_npl_mng      ! Region ID of north pole management
  integer, public              :: ADM_rgnid_spl_mng      ! Region ID of south pole management

  !
  !====== Information for regions ======
  !
  integer, public              :: ADM_l_me        ! Present Local region number

  logical, public, allocatable :: ADM_have_sgp(:) ! region have singlar point?

  !
  !====== Information for grids ======
  !
  integer, public, parameter   :: ADM_nxyz = 3 ! dimension of the spacial vector

  ! Identifiers of grid points around a grid point
  integer, public, parameter   :: ADM_GIJ_nmax = 7
  integer, public, parameter   :: ADM_GIoJo    = 1
  integer, public, parameter   :: ADM_GIpJo    = 2
  integer, public, parameter   :: ADM_GIpJp    = 3
  integer, public, parameter   :: ADM_GIoJp    = 4
  integer, public, parameter   :: ADM_GImJo    = 5
  integer, public, parameter   :: ADM_GImJm    = 6
  integer, public, parameter   :: ADM_GIoJm    = 7

#ifdef _FIXEDINDEX_
  integer, public              :: ADM_IooJoo(ADM_IooJoo_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_IooJmo(ADM_IooJmo_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_IooJop(ADM_IooJop_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_IooJmp(ADM_IooJmp_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_ImoJoo(ADM_ImoJoo_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_ImoJmo(ADM_ImoJmo_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_ImoJop(ADM_ImoJop_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_ImoJmp(ADM_ImoJmp_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_IopJoo(ADM_IopJoo_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_IopJmo(ADM_IopJmo_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_IopJop(ADM_IopJop_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_IopJmp(ADM_IopJmp_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_ImpJoo(ADM_ImpJoo_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_ImpJmo(ADM_ImpJmo_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_ImpJop(ADM_ImpJop_nmax,ADM_GIJ_nmax)
  integer, public              :: ADM_ImpJmp(ADM_ImpJmp_nmax,ADM_GIJ_nmax)
#else
  integer, public, allocatable :: ADM_IooJoo(:,:)
  integer, public, allocatable :: ADM_IooJmo(:,:)
  integer, public, allocatable :: ADM_IooJop(:,:)
  integer, public, allocatable :: ADM_IooJmp(:,:)
  integer, public, allocatable :: ADM_ImoJoo(:,:)
  integer, public, allocatable :: ADM_ImoJmo(:,:)
  integer, public, allocatable :: ADM_ImoJop(:,:)
  integer, public, allocatable :: ADM_ImoJmp(:,:)
  integer, public, allocatable :: ADM_IopJoo(:,:)
  integer, public, allocatable :: ADM_IopJmo(:,:)
  integer, public, allocatable :: ADM_IopJop(:,:)
  integer, public, allocatable :: ADM_IopJmp(:,:)
  integer, public, allocatable :: ADM_ImpJoo(:,:)
  integer, public, allocatable :: ADM_ImpJmo(:,:)
  integer, public, allocatable :: ADM_ImpJop(:,:)
  integer, public, allocatable :: ADM_ImpJmp(:,:)
#endif

  character(len=ADM_MAXFNAME), public :: ADM_HGRID_SYSTEM = 'ICO' ! [XTMS] Horizontal Grid type
                                                          ! 'ICO'      icosahedral
                                                          ! 'ICO-XTMS' icosahedral but XTMS is used in oprt
                                                          ! 'LCP'      Lambert-cornial (including PSP)
                                                          ! 'MLCP'     Mercator+Lambert-cornial
                                                          ! 'MLCP-OLD' OLD vergion (only for s=1)

  integer,                     public :: ADM_XTMS_MLCP_S  = 1 ! [XTMS] Number of segment for MLCP

  logical, public :: ADM_debug = .false.

  !-----------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: input_mnginfo
  private :: output_info
  private :: setup_vtab

  !-----------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private :: ADM_run_type ! Run type (single or multi processes)

  !-----------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------
  !> MPI initialization
  subroutine ADM_proc_init( &
       rtype )
#ifdef JCUP
    use jsp_nicam, only: &
       jsp_n_init, &
       jsp_n_get_my_mpi
#endif
    implicit none

    integer, intent(in) :: rtype

#ifdef JCUP
    integer :: my_comm, my_group
#endif
    integer :: my_rank, prc_all
    integer :: ierr
    !---------------------------------------------------------------------

    ADM_run_type = rtype

    if ( rtype == ADM_MULTI_PRC ) then
#ifdef JCUP
       call jsp_n_init("./nhm_driver.cnf")
       call jsp_n_get_my_mpi(my_comm, my_group, prc_all, my_rank)
       ADM_COMM_WORLD = my_comm
#else
       call MPI_Init(ierr)
       call MPI_Comm_size(MPI_COMM_WORLD, prc_all, ierr)
       call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
       ADM_COMM_WORLD = MPI_COMM_WORLD
#endif
       ADM_mpi_alive  = .true.
    else
       prc_all = 1
       my_rank = 0
    endif

    ADM_prc_me = my_rank + 1
    ADM_prc_pl = 1

#ifdef _FIXEDINDEX_
    if ( ADM_prc_all /= prc_all ) then
       write(*,*) 'xxx Fixed prc_all is not match (fixed,requested): ', ADM_prc_all, prc_all
       stop
    endif
#else
    ADM_prc_all = prc_all
#endif

    return
  end subroutine ADM_proc_init

  !-----------------------------------------------------------------------
  !> Abort MPI process
  subroutine ADM_proc_stop
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------

    ! flush 1kbyte
    write(ADM_LOG_FID,'(32A32)') '                                '

    write(ADM_LOG_FID,*) '+++ Abort MPI'
    if ( ADM_prc_me == ADM_prc_run_master ) then
       write(*,*) '+++ Abort MPI'
    endif

    close(ADM_LOG_FID)
    close(ADM_CTL_FID)

    ! Abort MPI
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

    stop
  end subroutine ADM_proc_stop

  !-----------------------------------------------------------------------
  !> Finish MPI process
  subroutine ADM_proc_finish
#ifdef JCUP
    use jsp_nicam, only: &
       jsp_n_finish,         &
       jsp_n_is_io_coupled,  &
       jsp_n_is_coco_coupled
#endif
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------

    if ( ADM_run_type == ADM_MULTI_PRC ) then

       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) '+++ finalize MPI'
       call MPI_Barrier(ADM_COMM_WORLD,ierr)

#ifdef JCUP
       if (      jsp_n_is_io_coupled()   &
            .OR. jsp_n_is_coco_coupled() ) then
          call jsp_n_finish()
       else
          call MPI_Finalize(ierr)
       endif
#else
       call MPI_Finalize(ierr)
#endif

       write(ADM_LOG_FID,*) '*** MPI is peacefully finalized'
    else
       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) '+++ stop serial process.'
    endif

    close(ADM_LOG_FID)
    close(ADM_CTL_FID)

    return
  end subroutine ADM_proc_finish

  !-----------------------------------------------------------------------
  !> Setup
  subroutine ADM_setup( &
       param_fname, &
       msg_base     )
    use mod_misc, only: &
       MISC_make_idstr, &
       MISC_get_available_fid, &
       MISC_get_fid
    implicit none

    character(LEN=*), intent(in) :: param_fname ! namelist file name

    character(len=*), intent(in), optional :: msg_base ! output file for msg.pexxxxx file

    integer                     :: glevel      = -1
    integer                     :: rlevel      = -1
    integer                     :: vlayer      =  1
    character(LEN=ADM_MAXFNAME) :: rgnmngfname = ''

    namelist / ADMPARAM / &
        glevel,           & !--- grid division level
        rlevel,           & !--- region division level
        vlayer,           & !--- number of inner vertical layer
        rgnmngfname,      & !--- region management file name
        ADM_HGRID_SYSTEM, & !--- grid system (default: ico)  ! S.Iga100607
        ADM_vlink_nmax,   & !--- num of lines at PL          ! S.Iga100607
        ADM_XTMS_MLCP_S,  & !--- num of segment for MLCP     ! S.Iga100607
        ADM_debug

    integer :: nmax, dmd
    integer :: l, rgnid
    integer :: ierr

    character(LEN=ADM_MAXFNAME) :: fname
    character(LEN=ADM_MAXFNAME) :: msg
    !---------------------------------------------------------------------

    msg = 'msg'
    if( present(msg_base) ) msg = msg_base ! [add] H.Yashiro 20110701

    !--- open message file
    ADM_LOG_FID = MISC_get_available_fid()
    call MISC_make_idstr(fname,trim(msg),'pe',ADM_prc_me)
    open( unit = ADM_LOG_FID, &
          file = trim(fname), &
          form = 'formatted'  )

    write(ADM_LOG_FID,*) '############################################################'
    write(ADM_LOG_FID,*) '#                                                          #'
    write(ADM_LOG_FID,*) '#   NICAM : Nonhydrostatic ICosahedal Atmospheric Model    #'
    write(ADM_LOG_FID,*) '#                                                          #'
    write(ADM_LOG_FID,*) '############################################################'

    ADM_CTL_FID = MISC_get_fid( param_fname )
    !--- open control file
    open( unit   = ADM_CTL_FID,       &
          file   = trim(param_fname), &
          form   = 'formatted',       &
          status = 'old',             &
          iostat = ierr               )

    if ( ierr /= 0 ) then
       write(*,*) 'xxx Cannot open parameter control file!'
       write(*,*) 'xxx filename:', trim(param_fname)
       call ADM_proc_stop
    endif

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[adm]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=ADMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(*,          *) 'xxx Not found namelist! STOP.'
       write(ADM_LOG_FID,*) 'xxx Not found namelist! STOP.'
       call ADM_proc_stop
    elseif ( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=ADMPARAM)

    ! Error if glevel & rlevel are not defined
    if ( glevel < 1 ) then
       write(*          ,*) 'xxx glevel is not appropriate :', glevel
       write(ADM_LOG_FID,*) 'xxx glevel is not appropriate :', glevel
       call ADM_proc_stop
    endif
    if ( rlevel < 0 ) then
       write(*          ,*) 'xxx rlevel is not appropriate :', rlevel
       write(ADM_LOG_FID,*) 'xxx rlevel is not appropriate :', rlevel
       call ADM_proc_stop
    endif

    ADM_rgnmngfname = trim(rgnmngfname)

    if ( ADM_HGRID_SYSTEM == 'ICO' ) then
       ADM_vlink_nmax = 5
       dmd            = 10
    elseif( ADM_HGRID_SYSTEM == 'LCP' ) then
       if( ADM_vlink_nmax == -1 ) ADM_vlink_nmax = 6
       dmd            = 4 * ADM_vlink_nmax
    elseif( ADM_HGRID_SYSTEM == 'MLCP-OLD' ) then
       if( ADM_vlink_nmax == -1 ) ADM_vlink_nmax = 6
       dmd            = 2 * ADM_vlink_nmax
    elseif( ADM_HGRID_SYSTEM == 'MLCP' ) then
       if( ADM_vlink_nmax == -1 ) ADM_vlink_nmax = 6
       dmd            = (1+ADM_XTMS_MLCP_S)  * ADM_vlink_nmax
    elseif( ADM_HGRID_SYSTEM == 'PERIODIC-1DMD' ) then ! T.Ohno 110721
       ADM_vlink_nmax = 5
       dmd            = 1
       ADM_prc_pl     = -999
    elseif( ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' ) then ! M.Hara 110721
       ADM_vlink_nmax = 5
       dmd            = 1
       ADM_prc_pl     = -999
    elseif( ADM_HGRID_SYSTEM == 'ICO-XTMS' ) then
       ADM_vlink_nmax = 5
       dmd            = 10
    else
       write(*          ,*) 'xxx Name of ADM_HGRID_SYSTEM is wrong. STOP.'
       write(ADM_LOG_FID,*) 'xxx Name of ADM_HGRID_SYSTEM is wrong. STOP.'
       call ADM_proc_stop
    endif

#ifdef _FIXEDINDEX_
    if ( ADM_vlink_nmax /= 5 ) then
       write(*          ,*) 'xxx Sorry, fixed index is not implemented for XTMS. STOP.'
       write(ADM_LOG_FID,*) 'xxx Sorry, fixed index is not implemented for XTMS. STOP.'
       call ADM_proc_stop
    endif
#else
    ADM_gall_pl = ADM_vlink_nmax + 1
    ADM_gmax_pl = ADM_vlink_nmax + 1
#endif



#ifdef _FIXEDINDEX_
    if ( ADM_glevel /= glevel ) then
       write(*,          *) 'xxx Fixed glevel is not match (fixed,requested): ', ADM_glevel, glevel
       write(ADM_LOG_FID,*) 'xxx Fixed glevel is not match (fixed,requested): ', ADM_glevel, glevel
       call ADM_proc_stop
    endif
    if ( ADM_rlevel /= rlevel ) then
       write(*,          *) 'xxx Fixed rlevel is not match (fixed,requested): ', ADM_rlevel, rlevel
       write(ADM_LOG_FID,*) 'xxx Fixed rlevel is not match (fixed,requested): ', ADM_rlevel, rlevel
       call ADM_proc_stop
    endif
    if ( ADM_vlayer /= vlayer ) then
       write(*,          *) 'xxx Fixed vlayer is not match (fixed,requested): ', ADM_vlayer, vlayer
       write(ADM_LOG_FID,*) 'xxx Fixed vlayer is not match (fixed,requested): ', ADM_vlayer, vlayer
       call ADM_proc_stop
    endif
    if ( ADM_DMD /= dmd ) then
       write(*,          *) 'xxx Fixed dmd is not match (fixed,requested): ', ADM_DMD, dmd
       write(ADM_LOG_FID,*) 'xxx Fixed dmd is not match (fixed,requested): ', ADM_DMD, dmd
       call ADM_proc_stop
    endif
#else
    ADM_glevel   = glevel
    ADM_rlevel   = rlevel
    ADM_vlayer   = vlayer
    ADM_DMD      = dmd

    ADM_rgn_nmax = 2**ADM_rlevel * 2**ADM_rlevel * ADM_DMD
    ADM_lall     = ADM_rgn_nmax / ADM_prc_all

    nmax         = 2**(ADM_glevel-ADM_rlevel)
    ADM_gall_1d  = 1 + nmax + 1
    ADM_gmin     = 1 + 1
    ADM_gmax     = 1 + nmax

    ADM_gall     = ( 1+nmax+1 ) * ( 1+nmax+1 )
    ADM_gall_in  = (   nmax+1 ) * (   nmax+1 )

    if ( ADM_vlayer == 1 ) then
       ADM_kall = 1
       ADM_kmin = 1
       ADM_kmax = 1
    else
       ADM_kall = 1 + ADM_vlayer + 1
       ADM_kmin = 1 + 1
       ADM_kmax = 1 + ADM_vlayer
    endif
#endif



    call input_mnginfo( ADM_rgnmngfname )

    ADM_prc_npl           = ADM_prc_pl
    ADM_prc_spl           = ADM_prc_pl
    ADM_prc_nspl(ADM_NPL) = ADM_prc_npl
    ADM_prc_nspl(ADM_SPL) = ADM_prc_spl

    allocate( ADM_have_sgp(ADM_lall) )
    ADM_have_sgp(:) = .false.

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)
       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          ADM_have_sgp(l) = .true.
       endif
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       ADM_have_pl = .true.
    else
       ADM_have_pl = .false.
    endif

    ! 2010.4.26 M.Satoh; 2010.5.11 M.Satoh
    ! ADM_l_me: this spans from 1 to ADM_lall, if effective.
    ! Otherwise, ADM_l_me = 0 should be set. see mod_history
    ADM_l_me = 0

    !--- make suffix for list-vector loop.
    call ADM_mk_suffix

    call output_info

    return
  end subroutine ADM_setup

  !-----------------------------------------------------------------------
  !> Read mnginfo file
  subroutine input_mnginfo( fname )
    use mod_misc,  only :&
       MISC_get_available_fid
    implicit none

    character(len=ADM_MAXFNAME), intent(in) :: fname

    integer :: num_of_rgn !--- number of region

    namelist / rgn_info / &
         num_of_rgn

    integer :: rgnid                    !--- region ID
    integer :: sw(ADM_RID:ADM_DIR) = -1 !--- south-west region info
    integer :: nw(ADM_RID:ADM_DIR) = -1 !--- nouth-west region info
    integer :: ne(ADM_RID:ADM_DIR) = -1 !--- nouth-east region info
    integer :: se(ADM_RID:ADM_DIR) = -1 !--- south-east region info

    namelist / rgn_link_info / &
         rgnid, &
         sw,    &
         nw,    &
         ne,    &
         se

    integer :: num_of_proc !--- number of run-processes

    namelist /proc_info/ &
         num_of_proc

    integer :: peid                         !--- process ID
    integer :: num_of_mng                   !--- number of regions be managed
    integer :: mng_rgnid(PRC_RGN_NMAX) = -1 !--- managed region ID

    namelist /rgn_mng_info/ &
         peid,       &
         num_of_mng, &
         mng_rgnid

    integer :: fid, ierr
    integer :: l, m, n
    !---------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[mnginfo]/Category[common share]'

    fid = MISC_get_available_fid()
    open( unit   = fid,         &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

    !=> [add] H.Yashiro 20120611
    ! ERROR if filename are not defined
    if ( ierr /= 0 ) then
       write(ADM_LOG_FID,*) 'xxx mnginfo file is not found! STOP. ', trim(fname)
       call ADM_proc_stop
    endif
    !<= [add] H.Yashiro 20120611

    read(fid,nml=rgn_info)
    if ( num_of_rgn /= ADM_rgn_nmax ) then
       write(ADM_LOG_FID,*) 'xxx No match for region number! STOP.'
       write(ADM_LOG_FID,*) 'xxx ADM_rgn_nmax= ',ADM_rgn_nmax,' num_of_rgn=',num_of_rgn
       call ADM_proc_stop
    endif

    allocate( ADM_rgn_etab( ADM_RID:ADM_DIR, &
                            ADM_SW:ADM_SE,   &
                            ADM_rgn_nmax     ) )

    do l = 1, ADM_rgn_nmax
       read(fid,nml=rgn_link_info)

       ADM_rgn_etab(ADM_RID:ADM_DIR,ADM_SW,rgnid) = sw(ADM_RID:ADM_DIR)
       ADM_rgn_etab(ADM_RID:ADM_DIR,ADM_NW,rgnid) = nw(ADM_RID:ADM_DIR)
       ADM_rgn_etab(ADM_RID:ADM_DIR,ADM_NE,rgnid) = ne(ADM_RID:ADM_DIR)
       ADM_rgn_etab(ADM_RID:ADM_DIR,ADM_SE,rgnid) = se(ADM_RID:ADM_DIR)
    enddo

    read(fid,nml=proc_info)
    if ( ADM_prc_all /= num_of_proc ) then
       write(ADM_LOG_FID,*) ' xxx No match for  process number! STOP.'
       write(ADM_LOG_FID,*) ' xxx ADM_prc_all= ',ADM_prc_all,' num_of_proc=',num_of_proc
       call ADM_proc_stop
    endif

    if ( ADM_prc_all /= num_of_proc ) then
       write(ADM_LOG_FID,*) 'Msg : Sub[ADM_input_mngtab]/Mod[admin]'
       write(ADM_LOG_FID,*) ' --- No match for process number!'
       call ADM_proc_stop
    endif

    allocate( ADM_prc_rnum(ADM_prc_all)              )
    allocate( ADM_prc_tab (PRC_RGN_NMAX,ADM_prc_all) )
    allocate( ADM_rgn2prc (ADM_rgn_nmax)             )
    ADM_prc_tab = -1 ! [Fix] 11/06/30  T.Seiki, fill undefined value

    do m = 1, ADM_prc_all
       read(fid,nml=rgn_mng_info)

       ADM_prc_rnum(m)      = num_of_mng
       ADM_prc_tab (:,peid) = mng_rgnid(:)
       do n = 1, num_of_mng
          ADM_rgn2prc(mng_rgnid(n)) = peid
       enddo
    enddo

    call setup_vtab

    close(fid)

    return
  end subroutine input_mnginfo

  !-----------------------------------------------------------------------
  subroutine setup_vtab
    implicit none

    integer :: nrid(ADM_vlink_nmax)
    integer :: nvid(ADM_vlink_nmax)
    integer :: vnum

    integer :: l, k, ll, v
    !---------------------------------------------------------------------

    allocate( ADM_rgn_vnum( ADM_W:ADM_S, &
                            ADM_rgn_nmax ) )

    allocate( ADM_rgn_vtab( ADM_RID:ADM_DIR,&
                            ADM_W:ADM_S,    &
                            ADM_rgn_nmax,   &
                            ADM_vlink_nmax  ) )

    allocate( ADM_rgn_vtab_pl( ADM_RID:ADM_DIR, &
                               ADM_rgn_nmax_pl, &
                               ADM_vlink_nmax   ) )

    do l = 1, ADM_rgn_nmax
       do k = ADM_W, ADM_S
          call set_vinfo(vnum,nrid,nvid,l,k)

          ADM_rgn_vnum(k,l)           = vnum
          ADM_rgn_vtab(ADM_RID,k,l,:) = nrid(:)
          ADM_rgn_vtab(ADM_DIR,k,l,:) = nvid(:)
       enddo
    enddo

    do l = 1, ADM_rgn_nmax
       if ( ADM_rgn_vnum(ADM_N,l) == ADM_vlink_nmax ) then
          ll = l
          exit
       endif
    enddo
    ADM_rgnid_npl_mng = ll

    do v = 1, ADM_vlink_nmax
       ADM_rgn_vtab_pl(ADM_RID,ADM_NPL,v) = ADM_rgn_vtab(ADM_RID,ADM_N,ll,v)
       ADM_rgn_vtab_pl(ADM_DIR,ADM_NPL,v) = ADM_rgn_vtab(ADM_DIR,ADM_N,ll,v)
    enddo

    do l = 1, ADM_rgn_nmax
       if ( ADM_rgn_vnum(ADM_S,l) == ADM_vlink_nmax ) then
          ll = l
          exit
       endif
    enddo
    ADM_rgnid_spl_mng = ll

    do v = 1, ADM_vlink_nmax
       ADM_rgn_vtab_pl(ADM_RID,ADM_SPL,v) = ADM_rgn_vtab(ADM_RID,ADM_S,ll,v)
       ADM_rgn_vtab_pl(ADM_DIR,ADM_SPL,v) = ADM_rgn_vtab(ADM_DIR,ADM_S,ll,v)
    enddo

    return
  end subroutine setup_vtab

  !-----------------------------------------------------------------------
  subroutine set_vinfo( vert_num, nrgnid, nvertid, rgnid, vertid )
    implicit none

    integer,intent(out) :: vert_num
    integer,intent(out) :: nrgnid (:)
    integer,intent(out) :: nvertid(:)
    integer,intent(in)  :: rgnid
    integer,intent(in)  :: vertid

    integer :: eid, rid
    integer :: eid_new, rid_new
    !---------------------------------------------------------------------

    vert_num = 0

    rid = rgnid
    eid = vertid
    select case(vertid)
    case(ADM_W)
       eid = ADM_SW
    case(ADM_N)
       eid = ADM_NW
    case(ADM_E)
       eid = ADM_NE
    case(ADM_S)
       eid = ADM_SE
    endselect

    nvertid(:) = -1
    nrgnid (:) = -1
    do
       rid_new = ADM_rgn_etab(ADM_RID,eid,rid)
       eid_new = ADM_rgn_etab(ADM_DIR,eid,rid) - 1

       if( eid_new == 0 ) eid_new = 4
       rid = rid_new
       eid = eid_new

       vert_num = vert_num + 1

       nrgnid (vert_num) = rid
       nvertid(vert_num) = eid

       if( rid == rgnid ) exit
    enddo

    return
  end subroutine set_vinfo

  !-----------------------------------------------------------------------
  subroutine ADM_mk_suffix
    implicit none

    integer :: gall_in
    integer :: i, j, n
    !---------------------------------------------------------------------

    gall_in = ADM_gmax - ADM_gmin + 1

#ifndef _FIXEDINDEX_
    ADM_IooJoo_nmax = ( gall_in   ) * ( gall_in   )
    ADM_IooJmo_nmax = ( gall_in   ) * ( gall_in+1 )
    ADM_IooJop_nmax = ( gall_in   ) * ( gall_in+1 )
    ADM_IooJmp_nmax = ( gall_in   ) * ( gall_in+2 )
    ADM_ImoJoo_nmax = ( gall_in+1 ) * ( gall_in   )
    ADM_ImoJmo_nmax = ( gall_in+1 ) * ( gall_in+1 )
    ADM_ImoJop_nmax = ( gall_in+1 ) * ( gall_in+1 )
    ADM_ImoJmp_nmax = ( gall_in+1 ) * ( gall_in+2 )
    ADM_IopJoo_nmax = ( gall_in+1 ) * ( gall_in   )
    ADM_IopJmo_nmax = ( gall_in+1 ) * ( gall_in+1 )
    ADM_IopJop_nmax = ( gall_in+1 ) * ( gall_in+1 )
    ADM_IopJmp_nmax = ( gall_in+1 ) * ( gall_in+2 )
    ADM_ImpJoo_nmax = ( gall_in+2 ) * ( gall_in   )
    ADM_ImpJmo_nmax = ( gall_in+2 ) * ( gall_in+1 )
    ADM_ImpJop_nmax = ( gall_in+2 ) * ( gall_in+1 )
    ADM_ImpJmp_nmax = ( gall_in+2 ) * ( gall_in+2 )

    allocate( ADM_IooJoo(ADM_IooJoo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IooJmo(ADM_IooJmo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IooJop(ADM_IooJop_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IooJmp(ADM_IooJmp_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImoJoo(ADM_ImoJoo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImoJmo(ADM_ImoJmo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImoJop(ADM_ImoJop_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImoJmp(ADM_ImoJmp_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IopJoo(ADM_IopJoo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IopJmo(ADM_IopJmo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IopJop(ADM_IopJop_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IopJmp(ADM_IopJmp_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImpJoo(ADM_ImpJoo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImpJmo(ADM_ImpJmo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImpJop(ADM_ImpJop_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImpJmp(ADM_ImpJmp_nmax,ADM_GIJ_nmax) )
#endif

    !--- ADM_IooJoo
    n = 1
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       ADM_IooJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJmo
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin,   ADM_gmax
       ADM_IooJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJop
    n = 1
    do j = ADM_gmin, ADM_gmax+1
    do i = ADM_gmin, ADM_gmax
       ADM_IooJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJmp
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin,   ADM_gmax
       ADM_IooJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJoo
    n = 1
    do j = ADM_gmin,   ADM_gmax
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJmo
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJop
    n = 1
    do j = ADM_gmin,   ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJmp
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJoo
    n = 1
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax+1
       ADM_IopJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJmo
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin,   ADM_gmax+1
       ADM_IopJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJop
    n = 1
    do j = ADM_gmin, ADM_gmax+1
    do i = ADM_gmin, ADM_gmax+1
       ADM_IopJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJmp
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin, ADM_gmax+1
       ADM_IopJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJoo
    n = 1
    do j = ADM_gmin,   ADM_gmax
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJmo
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJop
    n = 1
    do j = ADM_gmin,   ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJmp
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    return
  end subroutine ADM_mk_suffix

  !-----------------------------------------------------------------------
  subroutine output_info
    implicit none

    integer :: n, k, m
    integer :: rgnid
    !---------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,'(1x,A)'   ) '====== Process management info. ======'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Total number of process           : ', ADM_prc_all
    write(ADM_LOG_FID,'(1x,A,I7)') '--- My Process rank                   : ', ADM_prc_me
    write(ADM_LOG_FID,'(1x,A)'   ) '====== Region/Grid topology info. ======'
    write(ADM_LOG_FID,'(1x,A,A)' ) '--- Grid sysytem                      : ', ADM_HGRID_SYSTEM
    write(ADM_LOG_FID,'(1x,A,I7)') '--- #  of diamond                     : ', ADM_DMD
    write(ADM_LOG_FID,'(1x,A)'   ) '====== Region management info. ======'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Region level (RL)                 : ', ADM_rlevel
    write(ADM_LOG_FID,'(1x,A,I7,3(A,I4),A)') '--- Total number of region            : ', ADM_rgn_nmax, &
                                             ' (', 2**ADM_rlevel, ' x', 2**ADM_rlevel, ' x', ADM_DMD, ' )'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- #  of region per process          : ', ADM_lall
    write(ADM_LOG_FID,'(1x,A)'   ) '--- ID of region in my process        : '
    write(ADM_LOG_FID,*) ADM_prc_tab(1:ADM_lall, ADM_prc_me)

    write(ADM_LOG_FID,'(1x,A,I7)') '--- Region ID, contains north pole    : ', ADM_rgnid_npl_mng
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Region ID, contains south pole    : ', ADM_rgnid_spl_mng
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Process rank, managing north pole : ', ADM_prc_npl
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Process rank, managing south pole : ', ADM_prc_spl
    write(ADM_LOG_FID,'(1x,A)'   ) '====== Grid management info. ======'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Grid level (GL)                   : ', ADM_glevel
    write(ADM_LOG_FID,'(1x,A,I7,2(A,I4),A,I7,A)') '--- Total number of grid (horizontal) : ',  &
                                                  4**(ADM_glevel-ADM_rlevel)*ADM_rgn_nmax, &
                                                  ' (', 2**(ADM_glevel-ADM_rlevel),         &
                                                  ' x', 2**(ADM_glevel-ADM_rlevel),         &
                                                  ' x', ADM_rgn_nmax, ' )'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Number of vertical layer          : ', ADM_kmax-ADM_kmin+1

    if ( ADM_debug ) then
       do n = 1, ADM_lall
          rgnid = ADM_prc_tab(n, ADM_prc_me)
          write(ADM_LOG_FID,*) ' --- Link information for region', rgnid

          write(ADM_LOG_FID,*) '     < edge link >   --- ( rgnid , edgid )'
          do k = ADM_SW, ADM_SE
             write(ADM_LOG_FID,*) '     (',rgnid,',',k,') -> ',         &
                                  '(', ADM_rgn_etab(ADM_RID,k,rgnid),   &
                                  ',', ADM_rgn_etab(ADM_DIR,k,rgnid), ')'
          enddo

          write(ADM_LOG_FID,*) '     < vertex link > --- ( rgnid , edgid )'
          do k = ADM_W, ADM_S
             write(ADM_LOG_FID,*) '     (',rgnid,',',k,') : ', ADM_rgn_vnum(k,rgnid), 'point link'
             do m = 1, ADM_rgn_vnum(k,rgnid)
                write(ADM_LOG_FID,*) '                -> ',                  &
                                     '(', ADM_rgn_vtab(ADM_RID,k,rgnid,m),   &
                                     ',', ADM_rgn_vtab(ADM_DIR,k,rgnid,m), ')'
             enddo
          enddo

       enddo

       write(ADM_LOG_FID,*) ' --- Table of corresponding between region ID and process ID'
       write(ADM_LOG_FID,*) '    region ID :  process ID'
       do n = 1, ADM_rgn_nmax
          write(ADM_LOG_FID,'(I13,I14)') n, ADM_rgn2prc(n)
       enddo
    endif

    return
  end subroutine output_info

  !-----------------------------------------------------------------------------
  !> Get MPI time
  !> @return time
  function ADM_MPItime() result(time)
    implicit none

    real(RP) :: time
    !---------------------------------------------------------------------------

    if ( ADM_mpi_alive ) then
       time = real(MPI_WTIME(),kind=RP)
    else
       call cpu_time(time)
    endif

  end function ADM_MPItime

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_adm
