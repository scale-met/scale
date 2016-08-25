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
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ADM_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !====== Basic definition & information ======

  ! Fist colomn on the table for region and direction
  integer,  public, parameter :: ADM_RID = 1
  integer,  public, parameter :: ADM_DIR = 2

  ! Identifiers of directions of region edges
  integer,  public, parameter :: ADM_SW = 1
  integer,  public, parameter :: ADM_NW = 2
  integer,  public, parameter :: ADM_NE = 3
  integer,  public, parameter :: ADM_SE = 4

  ! Identifiers of directions of region vertices
  integer,  public, parameter :: ADM_W = 1
  integer,  public, parameter :: ADM_N = 2
  integer,  public, parameter :: ADM_E = 3
  integer,  public, parameter :: ADM_S = 4

  ! Identifier of triangle element (i-axis-side or j-axis side)
  integer,  public, parameter :: ADM_TI = 1
  integer,  public, parameter :: ADM_TJ = 2

  ! Identifier of arc element (i-axis-side, ij-axis side, or j-axis side)
  integer,  public, parameter :: ADM_AI  = 1
  integer,  public, parameter :: ADM_AIJ = 2
  integer,  public, parameter :: ADM_AJ  = 3

  ! Identifier of 1 variable
  integer,  public, parameter :: ADM_KNONE = 1

  ! Identifier of poles (north pole or south pole)
  integer,  public, parameter :: ADM_NPL = 1
  integer,  public, parameter :: ADM_SPL = 2

  ! dimension of the spacial vector
  integer,  public, parameter :: ADM_nxyz = 3

#ifdef _FIXEDINDEX_
  include "inc_index.h"
#else
  !#############################################################################
  ! Basic Index Parameters
  !#############################################################################

  ! main parameter
  integer,  public            :: ADM_glevel           ! grid   division level
  integer,  public            :: ADM_rlevel           ! region division level
  integer,  public            :: ADM_vlayer           ! number of vertical layer
  integer,  public            :: ADM_DMD              ! number of diamond

  ! region
  integer,  public            :: ADM_rgn_nmax         ! number of regular region
  integer,  public            :: ADM_lall             ! number of regular region per process
  integer,  public, parameter :: ADM_rgn_nmax_pl =  2 ! number of pole    region
  integer,  public, parameter :: ADM_lall_pl     =  2 ! number of pole    region per process

  ! horizontal grid
  integer,  public            :: ADM_gall             ! number of horizontal grid per regular region
  integer,  public            :: ADM_gall_in          ! number of horizontal grid (inner part)
  integer,  public            :: ADM_gall_1d          ! number of horizontal grid (1D)
  integer,  public            :: ADM_gmin             ! start index of 1D horizontal grid
  integer,  public            :: ADM_gmax             ! end   index of 1D horizontal grid

  integer,  public            :: ADM_iall             ! number of horizontal grid per regular region (i-axis)
  integer,  public            :: ADM_imin             ! start index of 1D horizontal grid
  integer,  public            :: ADM_imax             ! end   index of 1D horizontal grid
  integer,  public            :: ADM_jall             ! number of horizontal grid per regular region (j-axis)
  integer,  public            :: ADM_jmin             ! start index of 1D horizontal grid
  integer,  public            :: ADM_jmax             ! end   index of 1D horizontal grid

  integer,  public            :: ADM_vlink       = -1 ! maximum number of vertex linkage, ICO:5, PSP:6, LCP, MLCP:k
  integer,  public            :: ADM_gall_pl          ! number of horizontal grid for pole region
  integer,  public, parameter :: ADM_gslf_pl     =  1 ! index for pole point
  integer,  public, parameter :: ADM_gmin_pl     =  2 ! start index of grid around the pole point
  integer,  public            :: ADM_gmax_pl          ! end   index of grid around the pole point

  ! vertical grid
  integer,  public            :: ADM_kall             ! number of vertical grid
  integer,  public            :: ADM_kmin             ! start index of vertical grid
  integer,  public            :: ADM_kmax             ! end   index of vertical grid
#endif

  !
  !====== Information for processes ======
  !

  integer,  public            :: ADM_prc_me               ! my process ID
  integer,  public, parameter :: ADM_prc_master = 1       ! master process ID
  integer,  public            :: ADM_prc_pl               ! process ID which manages the pole regions
  logical,  public            :: ADM_have_pl              ! this ID manages pole region?

  integer,  public            :: ADM_prc_npl              ! process ID which have the pole regions
  integer,  public            :: ADM_prc_spl              ! process ID which have the pole regions
  integer,  public            :: ADM_prc_nspl(ADM_NPL:ADM_SPL)

  !
  !====== Information for processes-region relationship ======
  !
  character(len=H_LONG), public :: ADM_rgnmngfname        ! file name for region management info

  integer,  public, parameter   :: ADM_l_limit = 2560     ! maximum number of region per process.

  integer,  public, allocatable :: ADM_prc_rnum(:)        ! number of regions managed by each process = ADM_lall
  integer,  public, allocatable :: ADM_prc_tab (:,:)      ! table  of regions managed by each process

  integer,  public, allocatable :: ADM_rgn2prc (:)        ! Table of process ID from region ID
                                                          ! ADM_rgn2prc(ADM_rgn_nmax)
  integer,  public, allocatable :: ADM_rgn_etab(:,:,:)    ! table  of edge link information
                                                          ! ADM_rgn_etab( ADM_RID:ADM_DIR,ADM_SW:ADM_SE,ADM_rgn_nmax )
  integer,  public, allocatable :: ADM_rgn_vnum(:,:)      ! Table of n-vertex-link at the region vertex
                                                          ! ADM_rgn_vnum( ADM_W:ADM_S,ADM_rgn_nmax )
  integer,  public, allocatable :: ADM_rgn_vtab(:,:,:,:)  ! Table of vertex link information
                                                          ! ADM_rgn_vtab   ( ADM_RID:ADM_DIR, &
                                                          !                  ADM_W:ADM_S,     &
                                                          !                  ADM_rgn_nmax,    &
                                                          !                  ADM_vlink        )
  integer,  public, allocatable :: ADM_rgn_vtab_pl(:,:,:) ! Table of vertex link information for poles
                                                          ! ADM_rgn_vtab_pl( ADM_RID:ADM_DIR, &
                                                          !                  ADM_rgn_nmax_pl, &
                                                          !                  ADM_vlink        )

  integer,  public              :: ADM_rgnid_npl_mng      ! Region ID of north pole management
  integer,  public              :: ADM_rgnid_spl_mng      ! Region ID of south pole management

  !
  !====== Information for regions ======
  !
  integer,  public              :: ADM_l_me               ! Present Local region number

  logical,  public, allocatable :: ADM_have_sgp(:)        ! region have singlar point?

  !
  !====== Information for grids ======
  !
  character(len=H_SHORT), public :: ADM_HGRID_SYSTEM = 'ICO' ! [XTMS] Horizontal Grid type
                                                     ! 'ICO'      icosahedral
                                                     ! 'ICO-XTMS' icosahedral but XTMS is used in oprt
                                                     ! 'LCP'      Lambert-cornial (including PSP)
                                                     ! 'MLCP'     Mercator+Lambert-cornial
                                                     ! 'MLCP-OLD' OLD vergion (only for s=1)

  integer,               public :: ADM_XTMS_MLCP_S  = 1 ! [XTMS] Number of segment for MLCP

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: input_mnginfo
  private :: output_info
  private :: setup_vtab

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: ADM_debug = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ADM_setup
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_myrank,  &
       PRC_nprocs
    implicit none

    integer               :: glevel      = -1 !> grid division level
    integer               :: rlevel      = -1 !> region division level
    integer               :: vlayer      =  1 !> number of inner vertical layer
    character(LEN=H_LONG) :: rgnmngfname = '' !> region management file name

    namelist / ADMPARAM / &
        glevel,           &
        rlevel,           &
        vlayer,           &
        rgnmngfname,      &
        ADM_HGRID_SYSTEM, &
#ifndef _FIXEDINDEX_
        ADM_vlink,        &
#endif
        ADM_XTMS_MLCP_S,  &
        ADM_debug

    integer :: nmax, dmd
    integer :: l, rgnid
    integer :: ierr
    !---------------------------------------------------------------------------

    ADM_prc_me = PRC_myrank + 1
    ADM_prc_pl = 1

#ifdef _FIXEDINDEX_
    if ( ADM_prc_all /= PRC_nprocs ) then
       write(*,*) 'xxx Fixed prc_all is not match (fixed,requested): ', ADM_prc_all, PRC_nprocs
       stop
    endif
#endif

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[adm]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=ADMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(*         ,*) 'xxx Not found namelist! STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not found namelist! STOP.'
       call PRC_MPIstop
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=ADMPARAM)

    ! Error if glevel & rlevel are not defined
    if ( glevel < 1 ) then
       write(*         ,*) 'xxx glevel is not appropriate :', glevel
       if( IO_L ) write(IO_FID_LOG,*) 'xxx glevel is not appropriate :', glevel
       call PRC_MPIstop
    endif
    if ( rlevel < 0 ) then
       write(*         ,*) 'xxx rlevel is not appropriate :', rlevel
       if( IO_L ) write(IO_FID_LOG,*) 'xxx rlevel is not appropriate :', rlevel
       call PRC_MPIstop
    endif

    ADM_rgnmngfname = trim(rgnmngfname)

#ifdef _FIXEDINDEX_
    if ( ADM_HGRID_SYSTEM == 'ICO' ) then
       dmd        = 10
    elseif( ADM_HGRID_SYSTEM == 'PERIODIC-1DMD' ) then ! T.Ohno 110721
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' ) then ! M.Hara 110721
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == 'ICO-XTMS' ) then
       dmd        = 10
    else
       write(*         ,*) 'xxx Not appropriate param for ADM_HGRID_SYSTEM. STOP.', trim(ADM_HGRID_SYSTEM)
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not appropriate param for ADM_HGRID_SYSTEM. STOP.', trim(ADM_HGRID_SYSTEM)
       call PRC_MPIstop
    endif
#else
    if ( ADM_HGRID_SYSTEM == 'ICO' ) then
       ADM_vlink  = 5
       dmd        = 10
    elseif( ADM_HGRID_SYSTEM == 'LCP' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = 4 * ADM_vlink
    elseif( ADM_HGRID_SYSTEM == 'MLCP-OLD' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = 2 * ADM_vlink
    elseif( ADM_HGRID_SYSTEM == 'MLCP' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = (1+ADM_XTMS_MLCP_S)  * ADM_vlink
    elseif( ADM_HGRID_SYSTEM == 'PERIODIC-1DMD' ) then ! T.Ohno 110721
       ADM_vlink  = 5
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' ) then ! M.Hara 110721
       ADM_vlink  = 5
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == 'ICO-XTMS' ) then
       ADM_vlink  = 5
       dmd        = 10
    else
       write(*         ,*) 'xxx Not appropriate param for ADM_HGRID_SYSTEM. STOP.', trim(ADM_HGRID_SYSTEM)
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not appropriate param for ADM_HGRID_SYSTEM. STOP.', trim(ADM_HGRID_SYSTEM)
       call PRC_MPIstop
    endif

    ADM_gall_pl = ADM_vlink + 1
    ADM_gmax_pl = ADM_vlink + 1
#endif

#ifdef _FIXEDINDEX_
    if ( ADM_glevel /= glevel ) then
       write(*,         *) 'xxx Fixed glevel is not match (fixed,requested): ', ADM_glevel, glevel
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Fixed glevel is not match (fixed,requested): ', ADM_glevel, glevel
       call PRC_MPIstop
    endif
    if ( ADM_rlevel /= rlevel ) then
       write(*         ,*) 'xxx Fixed rlevel is not match (fixed,requested): ', ADM_rlevel, rlevel
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Fixed rlevel is not match (fixed,requested): ', ADM_rlevel, rlevel
       call PRC_MPIstop
    endif
    if ( ADM_vlayer /= vlayer ) then
       write(*         ,*) 'xxx Fixed vlayer is not match (fixed,requested): ', ADM_vlayer, vlayer
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Fixed vlayer is not match (fixed,requested): ', ADM_vlayer, vlayer
       call PRC_MPIstop
    endif
    if ( ADM_DMD /= dmd ) then
       write(*         ,*) 'xxx Fixed dmd is not match (fixed,requested): ', ADM_DMD, dmd
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Fixed dmd is not match (fixed,requested): ', ADM_DMD, dmd
       call PRC_MPIstop
    endif
#else
    ADM_glevel   = glevel
    ADM_rlevel   = rlevel
    ADM_vlayer   = vlayer
    ADM_DMD      = dmd

    ADM_rgn_nmax = 2**ADM_rlevel * 2**ADM_rlevel * ADM_DMD
    ADM_lall     = ADM_rgn_nmax / PRC_nprocs

    nmax         = 2**(ADM_glevel-ADM_rlevel)
    ADM_gall_1d  = 1 + nmax + 1
    ADM_gmin     = 1 + 1
    ADM_gmax     = 1 + nmax

    ADM_gall     = ( 1+nmax+1 ) * ( 1+nmax+1 )
    ADM_gall_in  = (   nmax+1 ) * (   nmax+1 )

    ADM_imin     = 1 + 1
    ADM_imax     = 1 + nmax
    ADM_iall     = 1 + nmax + 1
    ADM_jmin     = 1 + 1
    ADM_jmax     = 1 + nmax
    ADM_jall     = 1 + nmax + 1

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

    ADM_l_me = 0

    call output_info

    return
  end subroutine ADM_setup

  !-----------------------------------------------------------------------------
  subroutine input_mnginfo( fname )
    use scale_process, only: &
       PRC_nprocs, &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: fname

    integer :: num_of_rgn ! number of region

    namelist / rgn_info / &
         num_of_rgn

    integer :: rgnid                    ! region ID
    integer :: sw(ADM_RID:ADM_DIR) = -1 ! south-west region info
    integer :: nw(ADM_RID:ADM_DIR) = -1 ! nouth-west region info
    integer :: ne(ADM_RID:ADM_DIR) = -1 ! nouth-east region info
    integer :: se(ADM_RID:ADM_DIR) = -1 ! south-east region info

    namelist / rgn_link_info / &
         rgnid, &
         sw,    &
         nw,    &
         ne,    &
         se

    integer :: num_of_proc ! number of run-processes

    namelist /proc_info/ &
         num_of_proc

    integer :: peid                         ! process ID
    integer :: num_of_mng                   ! number of regions be managed
    integer :: mng_rgnid(ADM_l_limit) = -1 ! managed region ID

    namelist /rgn_mng_info/ &
         peid,       &
         num_of_mng, &
         mng_rgnid

    integer :: fid, ierr
    integer :: l, m, n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[mnginfo]/Category[common share]'

    fid = IO_get_available_fid()
    open( unit   = fid,         &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

    !=> [add] H.Yashiro 20120611
    ! ERROR if filename are not defined
    if ( ierr /= 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx mnginfo file is not found! STOP. ', trim(fname)
       call PRC_MPIstop
    endif
    !<= [add] H.Yashiro 20120611

    read(fid,nml=rgn_info)
    if ( num_of_rgn /= ADM_rgn_nmax ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx No match for region number! STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx ADM_rgn_nmax= ',ADM_rgn_nmax,' num_of_rgn=',num_of_rgn
       call PRC_MPIstop
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
    if ( PRC_nprocs /= num_of_proc ) then
       if( IO_L ) write(IO_FID_LOG,*) ' xxx No match for  process number! STOP.'
       if( IO_L ) write(IO_FID_LOG,*) ' xxx PRC_nprocs= ',PRC_nprocs,' num_of_proc=',num_of_proc
       call PRC_MPIstop
    endif

    if ( PRC_nprocs /= num_of_proc ) then
       if( IO_L ) write(IO_FID_LOG,*) 'Msg : Sub[ADM_input_mngtab]/Mod[admin]'
       if( IO_L ) write(IO_FID_LOG,*) ' --- No match for process number!'
       call PRC_MPIstop
    endif

    allocate( ADM_prc_rnum(PRC_nprocs)              )
    allocate( ADM_prc_tab (ADM_l_limit,PRC_nprocs) )
    allocate( ADM_rgn2prc (ADM_rgn_nmax)            )
    ADM_prc_tab = -1

    do m = 1, PRC_nprocs
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

  !-----------------------------------------------------------------------------
  subroutine setup_vtab
    implicit none

    integer :: nrid(ADM_vlink)
    integer :: nvid(ADM_vlink)
    integer :: vnum

    integer :: l, k, ll, v
    !---------------------------------------------------------------------------

    allocate( ADM_rgn_vnum( ADM_W:ADM_S, &
                            ADM_rgn_nmax ) )

    allocate( ADM_rgn_vtab( ADM_RID:ADM_DIR,&
                            ADM_W:ADM_S,    &
                            ADM_rgn_nmax,   &
                            ADM_vlink       ) )

    allocate( ADM_rgn_vtab_pl( ADM_RID:ADM_DIR, &
                               ADM_rgn_nmax_pl, &
                               ADM_vlink        ) )

    do l = 1, ADM_rgn_nmax
       do k = ADM_W, ADM_S
          call set_vinfo(vnum,nrid,nvid,l,k)

          ADM_rgn_vnum(k,l)           = vnum
          ADM_rgn_vtab(ADM_RID,k,l,:) = nrid(:)
          ADM_rgn_vtab(ADM_DIR,k,l,:) = nvid(:)
       enddo
    enddo

    do l = 1, ADM_rgn_nmax
       if ( ADM_rgn_vnum(ADM_N,l) == ADM_vlink ) then
          ll = l
          exit
       endif
    enddo
    ADM_rgnid_npl_mng = ll

    do v = 1, ADM_vlink
       ADM_rgn_vtab_pl(ADM_RID,ADM_NPL,v) = ADM_rgn_vtab(ADM_RID,ADM_N,ll,v)
       ADM_rgn_vtab_pl(ADM_DIR,ADM_NPL,v) = ADM_rgn_vtab(ADM_DIR,ADM_N,ll,v)
    enddo

    do l = 1, ADM_rgn_nmax
       if ( ADM_rgn_vnum(ADM_S,l) == ADM_vlink ) then
          ll = l
          exit
       endif
    enddo
    ADM_rgnid_spl_mng = ll

    do v = 1, ADM_vlink
       ADM_rgn_vtab_pl(ADM_RID,ADM_SPL,v) = ADM_rgn_vtab(ADM_RID,ADM_S,ll,v)
       ADM_rgn_vtab_pl(ADM_DIR,ADM_SPL,v) = ADM_rgn_vtab(ADM_DIR,ADM_S,ll,v)
    enddo

    return
  end subroutine setup_vtab

  !-----------------------------------------------------------------------------
  subroutine set_vinfo( vert_num, nrgnid, nvertid, rgnid, vertid )
    implicit none

    integer,intent(out) :: vert_num
    integer,intent(out) :: nrgnid (:)
    integer,intent(out) :: nvertid(:)
    integer,intent(in)  :: rgnid
    integer,intent(in)  :: vertid

    integer :: eid, rid
    integer :: eid_new, rid_new
    !---------------------------------------------------------------------------

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

  !-----------------------------------------------------------------------------
  subroutine output_info
    use scale_process, only: &
       PRC_nprocs
    implicit none

    integer :: n, k, m
    integer :: rgnid
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '====== Process management info. ======'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Total number of process           : ', PRC_nprocs
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- My Process rank                   : ', ADM_prc_me
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '====== Region/Grid topology info. ======'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)' ) '--- Grid sysytem                      : ', ADM_HGRID_SYSTEM
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- #  of diamond                     : ', ADM_DMD
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '====== Region management info. ======'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Region level (RL)                 : ', ADM_rlevel
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7,3(A,I4),A)') '--- Total number of region            : ', ADM_rgn_nmax, &
                                             ' (', 2**ADM_rlevel, ' x', 2**ADM_rlevel, ' x', ADM_DMD, ' )'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- #  of region per process          : ', ADM_lall
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '--- ID of region in my process        : '
    if( IO_L ) write(IO_FID_LOG,*)           ADM_prc_tab(1:ADM_lall, ADM_prc_me)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Region ID, contains north pole    : ', ADM_rgnid_npl_mng
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Region ID, contains south pole    : ', ADM_rgnid_spl_mng
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Process rank, managing north pole : ', ADM_prc_npl
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Process rank, managing south pole : ', ADM_prc_spl
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '====== Grid management info. ======'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Grid level (GL)                   : ', ADM_glevel
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7,2(A,I4),A,I7,A)') '--- Total number of grid (horizontal) : ', &
                                                 4**(ADM_glevel-ADM_rlevel)*ADM_rgn_nmax,    &
                                                 ' (', 2**(ADM_glevel-ADM_rlevel),           &
                                                 ' x', 2**(ADM_glevel-ADM_rlevel),           &
                                                 ' x', ADM_rgn_nmax, ' )'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Number of vertical layer          : ', ADM_kmax-ADM_kmin+1

    if ( ADM_debug ) then
       do n = 1, ADM_lall
          rgnid = ADM_prc_tab(n, ADM_prc_me)
          if( IO_L ) write(IO_FID_LOG,*) ' --- Link information for region', rgnid

          if( IO_L ) write(IO_FID_LOG,*) '     < edge link >   --- ( rgnid , edgid )'
          do k = ADM_SW, ADM_SE
             if( IO_L ) write(IO_FID_LOG,*) '     (',rgnid,',',k,') -> ',         &
                                 '(', ADM_rgn_etab(ADM_RID,k,rgnid),   &
                                 ',', ADM_rgn_etab(ADM_DIR,k,rgnid), ')'
          enddo

          if( IO_L ) write(IO_FID_LOG,*) '     < vertex link > --- ( rgnid , edgid )'
          do k = ADM_W, ADM_S
             if( IO_L ) write(IO_FID_LOG,*) '     (',rgnid,',',k,') : ', ADM_rgn_vnum(k,rgnid), 'point link'
             do m = 1, ADM_rgn_vnum(k,rgnid)
                if( IO_L ) write(IO_FID_LOG,*) '                -> ',                  &
                                    '(', ADM_rgn_vtab(ADM_RID,k,rgnid,m),   &
                                    ',', ADM_rgn_vtab(ADM_DIR,k,rgnid,m), ')'
             enddo
          enddo

       enddo

       if( IO_L ) write(IO_FID_LOG,*) ' --- Table of corresponding between region ID and process ID'
       if( IO_L ) write(IO_FID_LOG,*) '    region ID :  process ID'
       do n = 1, ADM_rgn_nmax
          if( IO_L ) write(IO_FID_LOG,'(I13,I14)') n, ADM_rgn2prc(n)
       enddo
    endif

    return
  end subroutine output_info

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_adm
