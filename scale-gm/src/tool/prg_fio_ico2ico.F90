!-------------------------------------------------------------------------------
!> Program FIO ico2ico
!!
!! @par Description
!!          This program converts between different resolution of icosahedral grids
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
program fio_ico2ico
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
     PRC_masterrank,       &
     PRC_MPIstart,         &
     PRC_SINGLECOM_setup,  &
     PRC_ERRHANDLER_setup, &
     PRC_MPIfinish,        &
     PRC_abort
  use scale_prc_icoA, only: &
     PRC_ICOA_RGN_generate, &
     I_l,                   &
     I_prc
  use scale_const, only: &
     UNDEF4 => CONST_UNDEF4, &
     UNDEF8 => CONST_UNDEF8
  use mod_io_param
  use mod_fio, only: &
       cstr, fstr
  use iso_c_binding
  !-----------------------------------------------------------------------------
  implicit none

  include 'fio_c.inc'
  !-----------------------------------------------------------------------------
  !
  !++ param & variable
  !
  integer, parameter :: max_nvar  = 500
  integer, parameter :: max_nstep = 1500

  !--- NAMELIST
  character(len=H_LONG)  :: src_fname           = ''
  integer                :: src_glevel          = -1
  integer                :: src_rlevel          = -1
!  integer                :: src_npe             = -1
  character(len=H_LONG)  :: dst_fname           = ''
  integer                :: dst_glevel          = -1
  integer                :: dst_rlevel          = -1
!  integer                :: dst_npe             = -1
  character(len=H_LONG)  :: iimap_fname         = ''

  logical                :: help = .false.

  namelist /PARAM_ICO2ICO/ &
     src_fname,      &
     src_PRC_RGN_ndiamond,   &
     src_glevel,     &
     src_rlevel,     &
     src_PRC_nprocs, &
     dst_fname,      &
     dst_PRC_RGN_ndiamond,   &
     dst_glevel,     &
     dst_rlevel,     &
     dst_PRC_nprocs, &
     iimap_fname,    &
     help

  !-----------------------------------------------------------------------------

  character(len=H_MID), parameter :: MODELNAME = "GM-ICO2ICO"
  character(len=H_LONG) :: cnf_fname
  integer               :: comm
  integer               :: nprocs
  integer               :: myrank
  logical               :: ismaster

  ! ico data information
  character(len=H_LONG) :: infname = ""
  integer, parameter    :: preclist(0:3) = (/ 4, 8, 4, 8 /)
  type(headerinfo) hinfo
  type(datainfo)   dinfo

  ! ico data
  integer               :: src_PRC_nprocs = -1
  integer               :: src_PRC_RGN_total
  integer               :: src_PRC_RGN_local
  integer               :: src_PRC_RGN_ndiamond = 10
  integer,  allocatable :: src_PRC_RGN_edge_tab(:,:,:) !< region link information (for 4 edges)
  integer,  allocatable :: src_PRC_RGN_lp2r    (:,:)   !< l,prc => rgn
  integer,  allocatable :: src_PRC_RGN_r2lp    (:,:)   !< rgn => l,prc
  integer               :: src_ADM_lall
  integer               :: src_ADM_gall
  integer               :: src_ADM_gall_1d
  integer               :: src_ADM_gmin
  integer               :: src_ADM_gmax
  integer,  allocatable :: src_prc_tab(:)
  integer,  allocatable :: src_prc_flag(:,:)
  integer,  allocatable :: src_fid     (:)
  real(SP), allocatable, target :: src_data4_1D(:)
  real(DP), allocatable, target :: src_data8_1D(:)
  real(SP), allocatable, target :: src_data4_3D(:,:,:)
  real(DP), allocatable, target :: src_data8_3D(:,:,:)
  integer               :: src_p, src_rgnid, src_g, src_l

  integer               :: dst_PRC_nprocs = -1
  integer               :: dst_PRC_RGN_total
  integer               :: dst_PRC_RGN_local
  integer               :: dst_PRC_RGN_ndiamond = 10  
  integer,  allocatable :: dst_PRC_RGN_edge_tab(:,:,:) !< region link information (for 4 edges)
  integer,  allocatable :: dst_PRC_RGN_lp2r    (:,:)   !< l,prc => rgn
  integer               :: dst_ADM_lall
  integer               :: dst_ADM_gall
  integer               :: dst_ADM_gall_1d
  integer               :: dst_ADM_gmin
  integer               :: dst_ADM_gmax
  integer,  allocatable :: dst_prc_tab(:)
  integer               :: dst_fid, dst_did
  real(SP), allocatable, target :: dst_data4_3D(:,:,:)
  real(DP), allocatable, target :: dst_data8_3D(:,:,:)
  integer               :: dst_p, dst_rgnid, dst_g, dst_l

  character(len=H_MID)  :: pkg_desc
  character(len=H_LONG) :: pkg_note
  integer               :: nmax_data
  integer               :: ADM_kall
  integer               :: did

  ! ii map
  real(DP), allocatable :: iimap(:,:,:,:,:)
  integer,  parameter   :: iimap_nmax = 7
  integer,  parameter   :: I_rgnid    = 1
  integer,  parameter   :: I_g1       = 2
  integer,  parameter   :: I_g2       = 3
  integer,  parameter   :: I_g3       = 4
  integer,  parameter   :: I_w1       = 5
  integer,  parameter   :: I_w2       = 6
  integer,  parameter   :: I_w3       = 7
  integer               :: g1, g2, g3
  real(DP)              :: weight
  real(DP), allocatable, target :: iibuf(:,:,:)


  ! for MPI
  integer :: prc_nlocal
  integer :: pstr, pend, pp

  character(len=IO_HSHORT) :: vname

  integer :: ierr
  integer :: k0 = 1
  integer :: k, firstfile
  !=============================================================================

  ! start MPI
  call PRC_MPIstart( comm ) ! [OUT]

  ! setup MPI communicator
  call PRC_SINGLECOM_setup( comm,    & ! [IN]
                            nprocs,  & ! [OUT]
                            myrank,  & ! [OUT]
                            ismaster ) ! [OUT]

  ! setup errhandler
  call PRC_ERRHANDLER_setup( .false., & ! [IN]
                             ismaster ) ! [IN]

  !########## Initial setup ##########

  ! setup standard I/O
  cnf_fname = IO_ARG_getfname( ismaster )

  call IO_setup( MODELNAME, cnf_fname )
  call IO_LOG_setup( myrank, ismaster )

  call PROF_rapstart('FIO_ICO2ICO_MPI')

  !--- read namelist
  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_ICO2ICO,iostat=ierr)
  if( ierr < 0 ) then !--- missing
     LOG_INFO("fio_ico2ico",*) 'Not found namelist. Default used.'
  elseif( ierr > 0 ) then !--- fatal error
     LOG_ERROR("fio_ico2ico",*) 'Not appropriate names in namelist PARAM_ICO2ICO. Check!'
     call PRC_abort
  endif
  LOG_NML(PARAM_ICO2ICO)

  if ( src_glevel==-1 .or. src_rlevel==-1 ) then
     LOG_ERROR("fio_ico2ico",*) "xxx Set src_glevel, src_rlevel. STOP"
     call PRC_abort
  endif
  if ( dst_glevel==-1 .or. dst_rlevel==-1 ) then
     LOG_ERROR("fio_ico2ico",*) "xxx Set dst_glevel, dst_rlevel. STOP"
     call PRC_abort
  endif
  if ( mod( dst_PRC_nprocs, nprocs ) /= 0 ) then
     LOG_ERROR("fio_ico2ico",*) "*** MPI processes must be divisible by dst_PRC_nprocs. STOP:", nprocs, dst_PRC_nprocs
     call PRC_abort
  endif

  !#########################################################

!  src_PRC_nprocs    = src_npe
  src_PRC_RGN_total = 10 * (4**src_rlevel)
  src_PRC_RGN_local = src_PRC_RGN_total / src_PRC_nprocs

  !--- prepare region infomation
  allocate( src_PRC_RGN_edge_tab(2,4,src_PRC_RGN_total) )
  allocate( src_PRC_RGN_lp2r    (src_PRC_RGN_local,0:src_PRC_nprocs-1) )
  allocate( src_PRC_RGN_r2lp    (2,src_PRC_RGN_total) )

  call PRC_ICOA_RGN_generate( src_rlevel,                  & ! [IN]
                              src_PRC_RGN_ndiamond,        & ! [IN] 
                              src_PRC_nprocs,              & ! [IN]
                              src_PRC_RGN_total,           & ! [IN]
                              src_PRC_RGN_local,           & ! [IN]
                              src_PRC_RGN_edge_tab(:,:,:), & ! [OUT]
                              src_PRC_RGN_lp2r    (:,:)    ) ! [OUT]

  do src_p = 0, src_PRC_nprocs-1
  do src_l = 1, src_PRC_RGN_local
     src_PRC_RGN_r2lp(I_l,  src_PRC_RGN_lp2r(src_l,src_p)) = src_l
     src_PRC_RGN_r2lp(I_prc,src_PRC_RGN_lp2r(src_l,src_p)) = src_p
  enddo
  enddo

  src_ADM_lall    = src_PRC_RGN_local
  src_ADM_gall_1d = 2**( src_glevel-src_rlevel ) + 2
  src_ADM_gmin    = 1               + 1
  src_ADM_gmax    = src_ADM_gall_1d - 1
  src_ADM_gall    = src_ADM_gall_1d * src_ADM_gall_1d

!  dst_PRC_nprocs    = dst_npe
  dst_PRC_RGN_total = 10 * (4**dst_rlevel)
  dst_PRC_RGN_local = dst_PRC_RGN_total / dst_PRC_nprocs

  !--- prepare region infomation
  allocate( dst_PRC_RGN_edge_tab(2,4,dst_PRC_RGN_total) )
  allocate( dst_PRC_RGN_lp2r    (dst_PRC_RGN_local,0:dst_PRC_nprocs-1) )

  call PRC_ICOA_RGN_generate( dst_rlevel,                  & ! [IN]
                              dst_PRC_RGN_ndiamond,        & ! [IN]
                              dst_PRC_nprocs,              & ! [IN]
                              dst_PRC_RGN_total,           & ! [IN]
                              dst_PRC_RGN_local,           & ! [IN]
                              dst_PRC_RGN_edge_tab(:,:,:), & ! [OUT]
                              dst_PRC_RGN_lp2r    (:,:)    ) ! [OUT]

  dst_ADM_lall    = dst_PRC_RGN_local
  dst_ADM_gall_1d = 2**( dst_glevel-dst_rlevel ) + 2
  dst_ADM_gmin    = 1               + 1
  dst_ADM_gmax    = dst_ADM_gall_1d - 1
  dst_ADM_gall    = dst_ADM_gall_1d * dst_ADM_gall_1d

  prc_nlocal = dst_PRC_nprocs / nprocs
  pstr       = myrank*prc_nlocal + 1
  pend       = myrank*prc_nlocal + prc_nlocal
  LOG_INFO("fio_ico2ico",*) "*** Number of Total .pexxxxxx files (source)     : ", src_PRC_nprocs
  LOG_INFO("fio_ico2ico",*) "*** Number of Total .pexxxxxx files (destination): ", dst_PRC_nprocs
  LOG_INFO("fio_ico2ico",*) "*** Number of PE to remapping process            : ", nprocs
  LOG_INFO("fio_ico2ico",*) "*** The rank of this process                     : ", myrank
  LOG_INFO("fio_ico2ico",*) "*** Number of files for this rank                : ", prc_nlocal
  LOG_INFO("fio_ico2ico",*) "*** file ID to pack                              : ", pstr-1, " - ", pend-1

  !--- setup
  ierr = fio_syscheck()

  !#########################################################

  LOG_INFO("fio_ico2ico",*) '*** iimap read start'
  call PROF_rapstart('READ IIMAP')

  allocate( iimap(dst_ADM_gall,k0,dst_ADM_lall,iimap_nmax,prc_nlocal) )
  allocate( iibuf(dst_ADM_gall,k0,dst_ADM_lall) )

  allocate( dst_prc_tab(dst_PRC_RGN_local) )

  do dst_p = pstr, pend
     pp = dst_p - pstr + 1

     dst_prc_tab(:) = dst_PRC_RGN_lp2r(:,dst_p-1)-1

     if ( pp == 1 ) then
        ierr = fio_put_commoninfo( IO_SPLIT_FILE,     &
                                   IO_BIG_ENDIAN,     &
                                   IO_ICOSAHEDRON,    &
                                   dst_glevel,        &
                                   dst_rlevel,        &
                                   dst_PRC_RGN_local, &
                                   dst_prc_tab        )
     endif

     call fio_mk_fname(infname,cstr(iimap_fname),cstr('pe'),dst_p-1,6)
     dst_fid = fio_register_file(cstr(infname))
     ierr = fio_fopen(dst_fid,IO_FREAD)
     ierr = fio_read_allinfo_validrgn(dst_fid,dst_prc_tab)

     did = fio_seek_datainfo(dst_fid,cstr("iimap_rgnid"),1)
     if ( did == -1 ) then
        write(*,*) 'xxx data not found! : iimap_rgnid'
        call PRC_abort
     endif
     ierr = fio_read_data(dst_fid,did,c_loc(iibuf))
     iimap(:,:,:,I_rgnid,pp) = iibuf(:,:,:)

     did = fio_seek_datainfo(dst_fid,cstr("iimap_g1"),1)
     ierr = fio_read_data(dst_fid,did,c_loc(iibuf))
     iimap(:,:,:,I_g1,pp) = iibuf(:,:,:)

     did = fio_seek_datainfo(dst_fid,cstr("iimap_g2"),1)
     ierr = fio_read_data(dst_fid,did,c_loc(iibuf))
     iimap(:,:,:,I_g2,pp) = iibuf(:,:,:)

     did = fio_seek_datainfo(dst_fid,cstr("iimap_g3"),1)
     ierr = fio_read_data(dst_fid,did,c_loc(iibuf))
     iimap(:,:,:,I_g3,pp) = iibuf(:,:,:)

     did = fio_seek_datainfo(dst_fid,cstr("iimap_w1"),1)
     ierr = fio_read_data(dst_fid,did,c_loc(iibuf))
     iimap(:,:,:,I_w1,pp) = iibuf(:,:,:)

     did = fio_seek_datainfo(dst_fid,cstr("iimap_w2"),1)
     ierr = fio_read_data(dst_fid,did,c_loc(iibuf))
     iimap(:,:,:,I_w2,pp) = iibuf(:,:,:)

     did = fio_seek_datainfo(dst_fid,cstr("iimap_w3"),1)
     ierr = fio_read_data(dst_fid,did,c_loc(iibuf))
     iimap(:,:,:,I_w3,pp) = iibuf(:,:,:)

     ierr = fio_fclose(dst_fid)
  enddo

  call PROF_rapend('READ IIMAP')
  LOG_INFO("fio_ico2ico",*) '*** iimap read end'

  !#########################################################

  LOG_INFO("fio_ico2ico",*) '*** icodata (source,header) read start'
  call PROF_rapstart('OPEN ICODATA')

  ! Check source file to read
  allocate( src_prc_flag(src_PRC_nprocs,prc_nlocal) )
  src_prc_flag(:,:) = 0

  do dst_p = pstr, pend
     pp = dst_p - pstr + 1

     do dst_l = 1, dst_ADM_lall
     do dst_g = 1, dst_ADM_gall
        dst_rgnid = iimap(dst_g,k0,dst_l,I_rgnid,pp)
        if ( dst_rgnid > 0 ) then
           src_p = src_PRC_RGN_r2lp(I_prc,dst_rgnid)

           src_prc_flag(src_p+1,pp) = 1
        endif
     enddo
     enddo
  enddo

  allocate( src_fid    (src_PRC_nprocs)    )
  allocate( src_prc_tab(src_PRC_RGN_local) )
  src_fid(:) = -1
  firstfile  = -1

  do src_p = 1, src_PRC_nprocs
     if( sum(src_prc_flag(src_p,:)) < 1 ) cycle

     src_prc_tab(:) = src_PRC_RGN_lp2r(:,src_p-1)-1

     if ( firstfile < 0 ) then
        ierr = fio_put_commoninfo( IO_SPLIT_FILE,     &
                                   IO_BIG_ENDIAN,     &
                                   IO_ICOSAHEDRON,    &
                                   src_glevel,        &
                                   src_rlevel,        &
                                   src_PRC_RGN_local, &
                                   src_prc_tab        )
     endif

     call fio_mk_fname(infname,cstr(src_fname),cstr('pe'),src_p-1,6)
     src_fid(src_p) = fio_register_file(cstr(infname))
     ierr = fio_fopen(src_fid(src_p),IO_FREAD)
     ierr = fio_read_allinfo_validrgn(src_fid(src_p),src_prc_tab)

     if ( firstfile < 0 ) then
        ierr = fio_get_pkginfo(hinfo,src_fid(src_p))
        call fstr(pkg_desc, hinfo%description)
        call fstr(pkg_note, hinfo%note)
        nmax_data = hinfo%num_of_data

        firstfile = src_p
     endif

  enddo

  call PROF_rapend('OPEN ICODATA')
  LOG_INFO("fio_ico2ico",*) '*** icodata (source,header) read end'

  !#########################################################

  LOG_INFO("fio_ico2ico",*) '*** convert start'
  call PROF_rapstart('CONVERT')

  do dst_p = pstr, pend
     pp = dst_p - pstr + 1

     dst_prc_tab(:) = dst_PRC_RGN_lp2r(:,dst_p-1)-1

     ierr = fio_put_commoninfo( IO_SPLIT_FILE,     &
                                IO_BIG_ENDIAN,     &
                                IO_ICOSAHEDRON,    &
                                dst_glevel,        &
                                dst_rlevel,        &
                                dst_PRC_RGN_local, &
                                dst_prc_tab        )

     call fio_mk_fname(infname,cstr(dst_fname),cstr('pe'),dst_p-1,6)
     dst_fid = fio_register_file(cstr(infname))
     ierr = fio_fopen(dst_fid,IO_FWRITE)
     ierr = fio_put_write_pkginfo(dst_fid,cstr(pkg_desc),cstr(pkg_note))

     LOG_INFO("fio_ico2ico",*) 'p(dst)=', dst_p, dst_fid

     do did = 0, nmax_data-1

        do src_p = 1, src_PRC_nprocs
           if( src_prc_flag(src_p,pp) < 1 ) cycle

           src_prc_tab(:) = src_PRC_RGN_lp2r(:,src_p-1)-1

           ierr = fio_put_commoninfo( IO_SPLIT_FILE,     &
                                      IO_BIG_ENDIAN,     &
                                      IO_ICOSAHEDRON,    &
                                      src_glevel,        &
                                      src_rlevel,        &
                                      src_PRC_RGN_local, &
                                      src_prc_tab        )

           if ( src_p == firstfile ) then
              ierr = fio_get_datainfo(dinfo,src_fid(src_p),did)
              ADM_kall = dinfo%num_of_layer
              call fstr(vname, dinfo%varname)

              if ( dinfo%datatype == IO_REAL4 ) then
                 allocate( src_data4_1D(src_ADM_gall*ADM_kall*src_ADM_lall) )
                 allocate( src_data4_3D(src_ADM_gall,ADM_kall,src_ADM_lall) )
                 allocate( dst_data4_3D(dst_ADM_gall,ADM_kall,dst_ADM_lall) )
                 dst_data4_3D(:,:,:) = 0.0_SP
              elseif( dinfo%datatype == IO_REAL8 ) then
                 allocate( src_data8_1D(src_ADM_gall*ADM_kall*src_ADM_lall) )
                 allocate( src_data8_3D(src_ADM_gall,ADM_kall,src_ADM_lall) )
                 allocate( dst_data8_3D(dst_ADM_gall,ADM_kall,dst_ADM_lall) )
                 dst_data8_3D(:,:,:) = 0.0_DP
              endif

              LOG_INFO("fio_ico2ico",*) '+varname:', dinfo%varname
              LOG_INFO("fio_ico2ico",*) ' +p(src)=', src_p, src_fid(src_p)
           endif

           if ( dinfo%datatype == IO_REAL4 ) then

              ierr = fio_read_data(src_fid(src_p),did,c_loc(src_data4_1D))
              src_data4_3D(:,:,:) = reshape( src_data4_1D(:), shape(src_data4_3D) )

              do dst_l = 1, dst_ADM_lall
              do dst_g = 1, dst_ADM_gall
                 dst_rgnid = nint( iimap(dst_g,k0,dst_l,I_rgnid,pp) )
                 g1        = nint( iimap(dst_g,k0,dst_l,I_g1   ,pp) )
                 g2        = nint( iimap(dst_g,k0,dst_l,I_g2   ,pp) )
                 g3        = nint( iimap(dst_g,k0,dst_l,I_g3   ,pp) )

                 do src_l = 1, src_ADM_lall
                    src_rgnid = src_PRC_RGN_lp2r(src_l,src_p-1)

                    if( src_rgnid /= dst_rgnid ) cycle

                    do src_g = 1, src_ADM_gall
                       if    ( src_g == g1 ) then
                          weight = iimap(dst_g,k0,dst_l,I_w1,pp)
                       elseif( src_g == g2 ) then
                          weight = iimap(dst_g,k0,dst_l,I_w2,pp)
                       elseif( src_g == g3 ) then
                          weight = iimap(dst_g,k0,dst_l,I_w3,pp)
                       else
                          weight = 0.0_DP
                       endif

                       if ( weight > 0.0_DP ) then
                          do k = 1, ADM_kall
                             if    ( src_data4_3D(src_g,k,src_l) <= UNDEF4 ) then
                                dst_data4_3D(dst_g,k,dst_l) = UNDEF4
                             else
                                dst_data4_3D(dst_g,k,dst_l) = dst_data4_3D(dst_g,k,dst_l) &
                                                            + src_data4_3D(src_g,k,src_l) * weight
                             endif
                          enddo
                       endif
                    enddo
                 enddo

              enddo
              enddo

           elseif( dinfo%datatype == IO_REAL8 ) then

              ierr = fio_read_data(src_fid(src_p),did,c_loc(src_data8_1D))
              src_data8_3D(:,:,:) = reshape( src_data8_1D(:), shape(src_data8_3D) )

              do dst_l = 1, dst_ADM_lall
              do dst_g = 1, dst_ADM_gall
                 dst_rgnid = nint( iimap(dst_g,k0,dst_l,I_rgnid,pp) )
                 g1        = nint( iimap(dst_g,k0,dst_l,I_g1   ,pp) )
                 g2        = nint( iimap(dst_g,k0,dst_l,I_g2   ,pp) )
                 g3        = nint( iimap(dst_g,k0,dst_l,I_g3   ,pp) )

                 do src_l = 1, src_ADM_lall
                    src_rgnid = src_PRC_RGN_lp2r(src_l,src_p-1)

                    if( src_rgnid /= dst_rgnid ) cycle

                    do src_g = 1, src_ADM_gall
                       if    ( src_g == g1 ) then
                          weight = iimap(dst_g,k0,dst_l,I_w1,pp)
                       elseif( src_g == g2 ) then
                          weight = iimap(dst_g,k0,dst_l,I_w2,pp)
                       elseif( src_g == g3 ) then
                          weight = iimap(dst_g,k0,dst_l,I_w3,pp)
                       else
                          weight = 0.0_DP
                       endif

                       if ( weight > 0.0_DP ) then
                          do k = 1, ADM_kall
                             if    ( src_data8_3D(src_g,k,src_l) <= UNDEF8 ) then
                                dst_data8_3D(dst_g,k,dst_l) = UNDEF8
                             else
                                dst_data8_3D(dst_g,k,dst_l) = dst_data8_3D(dst_g,k,dst_l) &
                                                         + src_data8_3D(src_g,k,src_l) * weight
                             endif
                          enddo
                       endif
                    enddo
                 enddo

              enddo
              enddo

           endif
        enddo

        ierr = fio_put_commoninfo( IO_SPLIT_FILE,     &
                                   IO_BIG_ENDIAN,     &
                                   IO_ICOSAHEDRON,    &
                                   dst_glevel,        &
                                   dst_rlevel,        &
                                   dst_PRC_RGN_local, &
                                   dst_prc_tab        )

        dinfo%datasize = int( dst_ADM_gall * dst_ADM_lall * ADM_kall * preclist(dinfo%datatype), kind=DP )

        if ( dinfo%datatype == IO_REAL4 ) then

           dst_did = fio_put_write_datainfo_data(dst_fid,dinfo,c_loc(dst_data4_3D))

           deallocate( src_data4_1D )
           deallocate( src_data4_3D )
           deallocate( dst_data4_3D )

        elseif( dinfo%datatype == IO_REAL8 ) then

           dst_did = fio_put_write_datainfo_data(dst_fid,dinfo,c_loc(dst_data8_3D))

           deallocate( src_data8_1D )
           deallocate( src_data8_3D )
           deallocate( dst_data8_3D )

        endif

     enddo

     ierr = fio_fclose(dst_fid)

  enddo

  do src_p = 1, src_PRC_nprocs
     if( sum(src_prc_flag(src_p,:)) < 1 ) cycle

     ierr = fio_fclose(src_fid(src_p))
  enddo

  call free( hinfo%rgnid )

  call PROF_rapend('CONVERT')
  LOG_INFO("fio_ico2ico",*) '*** convert finished! '

  call PROF_rapend('FIO_ICO2ICO_MPI')
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

end program fio_ico2ico
