!-------------------------------------------------------------------------------
!> Module mkiimap
!!
!! @par Description
!!         This module contains the tools to convert between two different icosahedral grid
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_mkiimap
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_atmos_grid_icoA_index
  use mod_io_param
  use iso_c_binding
  use mod_fio, only: cstr

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKIIMAP

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MKIIMAP_calcmap
  private :: MKIIMAP_setico
  private :: MKIIMAP_input_src_hgrid
  private :: MKIIMAP_output_iimap

  include 'fio_c.inc'
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private :: src_PRC_nprocs       = 0
  integer,                private :: src_PRC_RGN_level    = 0
  integer,                private :: src_PRC_RGN_ndiamond = 10  
  integer,                private :: src_GRID_LEVEL       = 0
  character(len=H_SHORT), private :: src_hgrid_io_mode    = 'ADVANCED'
  character(len=H_LONG),  private :: src_hgrid_fname      = ''
  character(len=H_SHORT), private :: OUT_io_mode          = 'ADVANCED'
  character(len=H_LONG),  private :: OUT_BASENAME         = ''
  real(RP),               private :: polar_limit          = -999.0_RP  ! search all longitude if abs(lat) > polar_limit

  integer,  private              :: src_PRC_RGN_total
  integer,  private              :: src_PRC_RGN_local
  integer,  private, allocatable :: src_PRC_RGN_edge_tab(:,:,:) !< region link information (for 4 edges)
  integer,  private, allocatable :: src_PRC_RGN_lp2r    (:,:)   !< l,prc => rgn

  integer,  private              :: src_ADM_lall
  integer,  private              :: src_ADM_gall
  integer,  private              :: src_ADM_gall_1d
  integer,  private              :: src_ADM_gmin
  integer,  private              :: src_ADM_gmax

  real(RP), private, allocatable :: src_GRD_x  (:,:,:,:)
  real(RP), private, allocatable :: src_GRD_LAT(:,:)
  real(RP), private, allocatable :: src_GRD_LON(:,:)

  integer,  private, allocatable :: checkmap(:,:,:,:)
  real(DP), private, allocatable :: iimap   (:,:,:,:)
  integer,  private, parameter   :: iimap_nmax = 7
  integer,  private, parameter   :: I_rgnid    = 1
  integer,  private, parameter   :: I_g1       = 2
  integer,  private, parameter   :: I_g2       = 3
  integer,  private, parameter   :: I_g3       = 4
  integer,  private, parameter   :: I_w1       = 5
  integer,  private, parameter   :: I_w2       = 6
  integer,  private, parameter   :: I_w3       = 7

  logical,  private :: debug = .false.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKIIMAP
    use scale_prc, only: &
       PRC_IsMaster, &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_RGN_level
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none

    real(RP) :: polar_limit_deg = 89.0_RP

    namelist / PARAM_MKIIMAP / &
       src_PRC_nprocs,       &
       src_PRC_RGN_level,    &
       src_PRC_RGN_ndiamond, &
       src_GRID_LEVEL,       &
       src_hgrid_io_mode,    &
       src_hgrid_fname,      &
       OUT_io_mode,          &
       OUT_BASENAME,         &
       polar_limit_deg,      &
       debug

    real(RP) :: weightsum

    integer  :: ierr
    integer  :: i, j, g, k, l, rgnid
    !---------------------------------------------------------------------------

    k = ADM_KNONE

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[mkiimap]/Category[prep]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKIIMAP,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** PARAM_MKIIMAP is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKIIMAP. STOP.'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_MKIIMAP)

    polar_limit = abs(polar_limit_deg) * D2R

    !--- setup source grid (icosahedral)
    call MKIIMAP_setico

    allocate( checkmap(ADM_gall,k,ADM_lall,1) )
    checkmap(:,:,:,:) = 0

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Start counting grid.'

    ! count lat-lon number
    call MKIIMAP_calcmap( 'GET_NUM' )

    if ( debug ) then
       do l = 1, ADM_lall
       do j = ADM_gmin, ADM_gmax+1
       do i = ADM_gmin, ADM_gmax+1
          g = suf(i,j)

          if ( checkmap(g,k,l,1) /= 1 ) then
             if( IO_L ) write(IO_FID_LOG,*) 'missed! (g,l,checkmap)=', &
                                            g, l, checkmap(g,k,l,1)
             call PRC_abort
          endif
       enddo
       enddo
       enddo
    endif

    allocate( iimap(ADM_gall,k,ADM_lall,iimap_nmax) )
    iimap(:,:,:,:) = -1.0_RP

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Start calc relation map.'

    ! calc relation map
    call MKIIMAP_calcmap( 'SET_INDEX' )

    if ( debug ) then
       do l = 1, ADM_lall
       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          g = suf(i,j)

          weightsum = iimap(g,k,l,I_w1) &
                    + iimap(g,k,l,I_w2) &
                    + iimap(g,k,l,I_w3)

          if ( abs(weightsum-1.0_RP) > 1.E-15_RP ) then
             if( IO_L ) write(IO_FID_LOG,'(A,2I6,E30.20)') 'invalid weight! (g,l,area)=', &
                                                           g, l, weightsum
          endif
       enddo
       enddo
       enddo
    endif

    ! output relation map
    call MKIIMAP_output_iimap

    return
  end subroutine MKIIMAP

  !-----------------------------------------------------------------------------
  subroutine MKIIMAP_calcmap( what_is_done )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
!       PRC_RGN_ndiamond,    &
       PRC_RGN_l2r,         &
       PRC_RGN_rgn4pl,      &
       I_NPL,               &
       I_SPL
    use scale_const, only: &
       PI => CONST_PI
    use scale_vector, only: &
       VECTR_triangle, &
       VECTR_cross,    &
       VECTR_dot,      &
       VECTR_abs
    use mod_grd, only: &
       GRD_rscale, &
       GRD_x,      &
       GRD_LAT,    &
       GRD_LON
    use mod_gmtr, only: &
       GMTR_polygon_type
    implicit none

    character(len=*), intent(in) :: what_is_done

    real(RP), parameter :: o(3) = 0.0_RP
    real(RP) :: r0(3), r1(3), r2(3), r3(3)
    real(RP) :: nvec(3), v12xv10(3), v23xv20(3), v31xv30(3)
    real(RP) :: judge12, judge23, judge31
    real(RP) :: rf, rn, ip, len
    real(RP) :: v01(3)

    real(RP) :: lat, lon
    real(RP) :: lat1, lat2, lat3
    real(RP) :: lon1, lon2, lon3
    real(RP) :: latmin_l,latmax_l
    real(RP) :: lonmin_l,lonmax_l
    logical  :: near_pole

    real(RP) :: area_total, area1, area2, area3

    real(RP) :: eps_judge  = 1.E-15_RP ! marginal value for inner products
    real(RP) :: eps_latlon             ! marginal square near grid points (in radian)
    real(RP) :: eps_vertex = 1.E-10_RP ! marginal value for vartex

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1

    integer :: rgnid
    integer :: i, j, g, k, l, isrc, jsrc, lsrc, t, p
    !---------------------------------------------------------------------------

    k = ADM_KNONE
    eps_latlon = sqrt( 4.0_RP * PI / real(src_PRC_RGN_ndiamond*4**src_GRID_LEVEL,kind=RP) )

    do p = 0, src_PRC_nprocs-1

       call MKIIMAP_input_src_hgrid( p )

       do lsrc = 1, src_ADM_lall
          rgnid = src_PRC_RGN_lp2r(lsrc,p)

          do jsrc = src_ADM_gmin, src_ADM_gmax
          do isrc = src_ADM_gmin, src_ADM_gmax

             ij     = src_ADM_gall_1d*(jsrc-1)+isrc
             ip1j   = src_ADM_gall_1d*(jsrc-1)+isrc+1
             ip1jp1 = src_ADM_gall_1d*(jsrc  )+isrc+1
             ijp1   = src_ADM_gall_1d*(jsrc  )+isrc

             do t = ADM_TI, ADM_TJ

                if ( t == ADM_TI ) then
                   r1(:) = src_GRD_x(ij    ,k,lsrc,:)
                   r2(:) = src_GRD_x(ip1j  ,k,lsrc,:)
                   r3(:) = src_GRD_x(ip1jp1,k,lsrc,:)

                   lat1 = src_GRD_LAT(ij    ,lsrc)
                   lon1 = src_GRD_LON(ij    ,lsrc)
                   lat2 = src_GRD_LAT(ip1j  ,lsrc)
                   lon2 = src_GRD_LON(ip1j  ,lsrc)
                   lat3 = src_GRD_LAT(ip1jp1,lsrc)
                   lon3 = src_GRD_LON(ip1jp1,lsrc)
                else !--- ADM_TJ
                   r1(:) = src_GRD_x(ij    ,k,lsrc,:)
                   r2(:) = src_GRD_x(ip1jp1,k,lsrc,:)
                   r3(:) = src_GRD_x(ijp1  ,k,lsrc,:)

                   lat1 = src_GRD_LAT(ij    ,lsrc)
                   lon1 = src_GRD_LON(ij    ,lsrc)
                   lat2 = src_GRD_LAT(ip1jp1,lsrc)
                   lon2 = src_GRD_LON(ip1jp1,lsrc)
                   lat3 = src_GRD_LAT(ijp1  ,lsrc)
                   lon3 = src_GRD_LON(ijp1  ,lsrc)
                endif

                latmax_l = max(lat1,lat2,lat3) + eps_latlon
                latmin_l = min(lat1,lat2,lat3) - eps_latlon

                if( latmin_l >  polar_limit ) latmax_l =  PI
                if( latmax_l < -polar_limit ) latmin_l = -PI

                lonmax_l = max(lon1,lon2,lon3)
                lonmin_l = min(lon1,lon2,lon3)
                if ( lonmax_l-lonmin_l > PI ) then
                   if( lon1 < 0 ) lon1 = lon1 + 2.0_RP * PI
                   if( lon2 < 0 ) lon2 = lon2 + 2.0_RP * PI
                   if( lon3 < 0 ) lon3 = lon3 + 2.0_RP * PI

                   lonmax_l = max(lon1,lon2,lon3)
                   lonmin_l = min(lon1,lon2,lon3)
                endif
                lonmax_l = lonmax_l + eps_latlon
                lonmin_l = lonmin_l - eps_latlon

                do l = 1, ADM_lall
                do j = ADM_gmin, ADM_gmax+1
                do i = ADM_gmin, ADM_gmax+1
                   g = suf(i,j)

                   lat = GRD_LAT(g,l)
                   lon = GRD_LON(g,l)

                   if( lat > latmax_l ) cycle
                   if( lat < latmin_l ) cycle

                   near_pole = .false.
                   if( lat >  polar_limit ) near_pole = .true.
                   if( lat < -polar_limit ) near_pole = .true.

                   if ( .NOT. near_pole ) then
                      if ( .NOT. (      (       ( lon             <= lonmax_l ) &
                                          .AND. ( lon             >= lonmin_l ) ) &
                                   .OR. (       ( lon - 2.0_RP*PI <= lonmax_l ) &
                                          .AND. ( lon - 2.0_RP*PI >= lonmin_l ) ) &
                                   .OR. (       ( lon + 2.0_RP*PI <= lonmax_l ) &
                                          .AND. ( lon + 2.0_RP*PI >= lonmin_l ) ) ) ) then
                         cycle
                      endif
                   endif

                   !--- target latlon point on the sphere
                   r0(:) = GRD_x(g,k,l,:) / GRD_rscale

                   !--- remove the case inner product is negative
                   call VECTR_dot( ip, o(:), r1(:), o(:), r0(:) )
                   if( ip < 0.0_RP ) cycle
                   v01(:) = r1(:) - r0(:)

                   !--- normal vector
                   call VECTR_cross( nvec(:), r1(:), r2(:), r2(:), r3(:) )
                   call VECTR_abs( len, nvec(:) )

                   nvec(:) = nvec(:) / len

                   !------ distance from origin to a plane with r0.
                   call VECTR_dot( rf, o(:), nvec(:), o(:), r0(:) )
                   !------ distance from origin to a plane with r1 or (r2,r3).
                   call VECTR_dot( rn, o(:), nvec(:), o(:), r1(:) )

                   !------ mapping r0
                   r0(1) = r0(1) * (rn/rf)
                   r0(2) = r0(2) * (rn/rf)
                   r0(3) = r0(3) * (rn/rf)

                   !--- calculate vectors from triangler points
                   call VECTR_cross( v12xv10(:), r1(:), r2(:), r0(:), r1(:) )
                   call VECTR_cross( v23xv20(:), r2(:), r3(:), r0(:), r2(:) )
                   call VECTR_cross( v31xv30(:), r3(:), r1(:), r0(:), r3(:) )

                   call VECTR_dot( judge12, o(:), nvec(:), o(:), v12xv10(:) )
                   call VECTR_dot( judge23, o(:), nvec(:), o(:), v23xv20(:) )
                   call VECTR_dot( judge31, o(:), nvec(:), o(:), v31xv30(:) )

                   if (       judge12 < eps_judge &
                        .AND. judge23 < eps_judge &
                        .AND. judge31 < eps_judge ) then ! in the triangle

                      select case( trim(what_is_done) )
                      case( 'GET_NUM' )

                         checkmap(g,k,l,1) = 1

                      case('SET_INDEX')

                         if( iimap(g,k,l,I_g1) > 0.0_RP ) cycle

                         iimap(g,k,l,I_rgnid) = rgnid
                         if ( t == ADM_TI ) then
                            iimap(g,k,l,I_g1) = ij
                            iimap(g,k,l,I_g2) = ip1j
                            iimap(g,k,l,I_g3) = ip1jp1
                         else !--- ADM_TJ
                            iimap(g,k,l,I_g1) = ij
                            iimap(g,k,l,I_g2) = ip1jp1
                            iimap(g,k,l,I_g3) = ijp1
                         endif

                         area1 = VECTR_triangle( r0(:), r2(:), r3(:),          &
                                                 GMTR_polygon_type, GRD_rscale )
                         area2 = VECTR_triangle( r0(:), r3(:), r1(:),          &
                                                 GMTR_polygon_type, GRD_rscale )
                         area3 = VECTR_triangle( r0(:), r1(:), r2(:),          &
                                                 GMTR_polygon_type, GRD_rscale )

                         if (      area1 * 0.0_RP /= 0.0_RP &
                              .OR. area2 * 0.0_RP /= 0.0_RP &
                              .OR. area3 * 0.0_RP /= 0.0_RP ) then ! Nan?
                            write(*,*) 'Nan! (g,l)=', g,l
                            write(*,*) '(rgnid,g1,area1)=', rgnid,iimap(g,k,l,I_g1),area1
                            write(*,*) '(rgnid,g2,area2)=', rgnid,iimap(g,k,l,I_g2),area2
                            write(*,*) '(rgnid,g3,area3)=', rgnid,iimap(g,k,l,I_g3),area3
                            if( IO_L ) write(IO_FID_LOG,*) 'Nan! (g,l)=', g,l
                            if( IO_L ) write(IO_FID_LOG,*) '(rgnid,g1,area1)=', rgnid,iimap(g,k,l,I_g1),area1
                            if( IO_L ) write(IO_FID_LOG,*) '(rgnid,g2,area2)=', rgnid,iimap(g,k,l,I_g2),area2
                            if( IO_L ) write(IO_FID_LOG,*) '(rgnid,g3,area3)=', rgnid,iimap(g,k,l,I_g3),area3
                            call PRC_abort
                         endif

                         area_total = area1 + area2 + area3

                         iimap(g,k,l,I_w1) = area1 / area_total
                         iimap(g,k,l,I_w2) = area2 / area_total
                         iimap(g,k,l,I_w3) = area3 / area_total

                      endselect

                      cycle

                   elseif(       t == ADM_TI              &
                           .AND. abs(v01(1)) < eps_vertex &
                           .AND. abs(v01(2)) < eps_vertex &
                           .AND. abs(v01(3)) < eps_vertex ) then ! on the triangle vertex

                      select case( trim(what_is_done) )
                      case( 'GET_NUM' )

                         checkmap(g,k,l,1) = 1

                      case('SET_INDEX')

                         if( iimap(g,k,l,I_g1) > 0.0_RP ) cycle

                         iimap(g,k,l,I_rgnid) = rgnid
                         iimap(g,k,l,I_g1)    = ij
                         iimap(g,k,l,I_g2)    = ip1j
                         iimap(g,k,l,I_g3)    = ip1jp1
                         iimap(g,k,l,I_w1)    = 1.0_RP
                         iimap(g,k,l,I_w2)    = 0.0_RP
                         iimap(g,k,l,I_w3)    = 0.0_RP

                      endselect

                   endif

                enddo ! i LOOP
                enddo ! j LOOP
                enddo ! l LOOP

             enddo ! TI,TJ
          enddo ! isrc LOOP
          enddo ! jsrc LOOP
       enddo ! lsrc LOOP
    enddo ! p LOOP

    if( IO_L ) write(IO_FID_LOG,*) 'OK.'

    return
  end subroutine MKIIMAP_calcmap

  !-----------------------------------------------------------------------------
  subroutine MKIIMAP_setico
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_ICOA_RGN_generate, &
!       PRC_RGN_ndiamond,      &
       I_RGNID, I_DIR,        &
       I_SW, I_SE
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    !---------------------------------------------------------------------------

    if ( src_PRC_RGN_level < 0 ) then
       write(IO_FID_LOG,*) 'Please set src_PRC_RGN_level! STOP.'
       call PRC_abort
    endif

    if ( src_GRID_LEVEL < 0 ) then
       write(IO_FID_LOG,*) 'Please set src_GRID_LEVEL! STOP.'
       call PRC_abort
    endif

    ! setup proc and region information
    src_PRC_RGN_total = 2**src_PRC_RGN_level * 2**src_PRC_RGN_level * src_PRC_RGN_ndiamond
    src_PRC_RGN_local = src_PRC_RGN_total / src_PRC_nprocs

    allocate( src_PRC_RGN_edge_tab(I_RGNID:I_DIR,I_SW:I_SE,src_PRC_RGN_total) )
    allocate( src_PRC_RGN_lp2r    (src_PRC_RGN_local,0:src_PRC_nprocs-1) )

    call PRC_ICOA_RGN_generate( src_PRC_RGN_level,           & ! [IN]
                                src_PRC_RGN_ndiamond,        & ! [IN]         
                                src_PRC_nprocs,              & ! [IN]
                                src_PRC_RGN_total,           & ! [IN]
                                src_PRC_RGN_local,           & ! [IN]
                                src_PRC_RGN_edge_tab(:,:,:), & ! [OUT]
                                src_PRC_RGN_lp2r    (:,:)    ) ! [OUT]

    ! setup grid information
    src_ADM_lall    = src_PRC_RGN_local
    src_ADM_gall_1d = 2**( src_GRID_LEVEL-src_PRC_RGN_level ) + 2
    src_ADM_gmin    = 1               + 1
    src_ADM_gmax    = src_ADM_gall_1d - 1
    src_ADM_gall    = src_ADM_gall_1d * src_ADM_gall_1d

    allocate( src_GRD_x  (src_ADM_gall,ADM_KNONE,src_ADM_lall,ADM_nxyz) )
    allocate( src_GRD_LAT(src_ADM_gall,src_ADM_lall) )
    allocate( src_GRD_LON(src_ADM_gall,src_ADM_lall) )
    src_GRD_x  (:,:,:,:) = UNDEF
    src_GRD_LAT(:,:)   = UNDEF
    src_GRD_LON(:,:)   = UNDEF

    return
  end subroutine MKIIMAP_setico

  !-----------------------------------------------------------------------------
  subroutine MKIIMAP_input_src_hgrid( &
       prcid )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_RGN_level, &
       PRC_RGN_local, &
       PRC_RGN_l2r
    use scale_vector, only: &
       VECTR_xyz2latlon
    use mod_grd, only: &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR
    implicit none

    integer, intent(in) :: prcid

    character(len=H_LONG) :: infname = ""
    integer, allocatable  :: prc_tab(:)
    type(datainfo) dinfo

    real(SP), target :: var4(src_ADM_gall,ADM_KNONE,src_ADM_lall)
    real(DP), target :: var8(src_ADM_gall,ADM_KNONE,src_ADM_lall)

    integer :: did, fid
    integer :: ierr
    integer :: g, k, l
    !---------------------------------------------------------------------------

    k = ADM_KNONE

    if ( src_hgrid_io_mode == 'ADVANCED' ) then

       allocate( prc_tab(src_PRC_RGN_local) )
       prc_tab(:) = src_PRC_RGN_lp2r(:,prcid)-1

       call fio_mk_fname(infname,cstr(src_hgrid_fname),cstr('pe'),prcid,6)

       if ( prcid == 0 ) then
          ierr = fio_put_commoninfo( IO_SPLIT_FILE,      & ! [IN]
                                     IO_BIG_ENDIAN,      & ! [IN]
                                     IO_ICOSAHEDRON,     & ! [IN]
                                     src_GRID_LEVEL,     & ! [IN]
                                     src_PRC_RGN_level,  & ! [IN]
                                     src_PRC_RGN_local,  & ! [IN]
                                     prc_tab             ) ! [IN]
       endif

       fid = fio_register_file(cstr(infname))
       ierr = fio_fopen(fid,IO_FREAD)
       ierr = fio_read_allinfo_validrgn(fid,prc_tab)

       deallocate( prc_tab )

       did = fio_seek_datainfo(fid,cstr("grd_x_x"),1)
       if ( did == -1 ) then
          write(*,*) 'xxx data not found! : grd_x_x'
          call PRC_abort
       endif
       ierr = fio_get_datainfo(dinfo,fid,did)

       if ( dinfo%datatype == IO_REAL4 ) then
          ierr = fio_read_data(fid,did,c_loc(var4))
          src_GRD_x(:,k,:,GRD_XDIR) = real( var4(:,k,:), kind=RP )
       elseif( dinfo%datatype == IO_REAL8 ) then
          ierr = fio_read_data(fid,did,c_loc(var8))
          src_GRD_x(:,k,:,GRD_XDIR) = real( var8(:,k,:), kind=RP )
       endif

       did = fio_seek_datainfo(fid,cstr("grd_x_y"),1)
       if ( did == -1 ) then
          write(*,*) 'xxx data not found! : grd_x_y'
          call PRC_abort
       endif
       ierr = fio_get_datainfo(dinfo,fid,did)

       if ( dinfo%datatype == IO_REAL4 ) then
          ierr = fio_read_data(fid,did,c_loc(var4))
          src_GRD_x(:,k,:,GRD_YDIR) = real( var4(:,k,:), kind=RP )
       elseif( dinfo%datatype == IO_REAL8 ) then
          ierr = fio_read_data(fid,did,c_loc(var8))
          src_GRD_x(:,k,:,GRD_YDIR) = real( var8(:,k,:), kind=RP )
       endif

       did = fio_seek_datainfo(fid,cstr("grd_x_z"),1)
       if ( did == -1 ) then
          write(*,*) 'xxx data not found! : grd_x_z'
          call PRC_abort
       endif
       ierr = fio_get_datainfo(dinfo,fid,did)

       if ( dinfo%datatype == IO_REAL4 ) then
          ierr = fio_read_data(fid,did,c_loc(var4))
          src_GRD_x(:,k,:,GRD_ZDIR) = real( var4(:,k,:), kind=RP )
       elseif( dinfo%datatype == IO_REAL8 ) then
          ierr = fio_read_data(fid,did,c_loc(var8))
          src_GRD_x(:,k,:,GRD_ZDIR) = real( var8(:,k,:), kind=RP )
       endif

       do l = 1, src_ADM_lall
       do g = 1, src_ADM_gall
          call VECTR_xyz2latlon( src_GRD_x(g,k,l,GRD_XDIR), & ! [IN]
                                 src_GRD_x(g,k,l,GRD_YDIR), & ! [IN]
                                 src_GRD_x(g,k,l,GRD_ZDIR), & ! [IN]
                                 src_GRD_LAT(g,l),          & ! [OUT]
                                 src_GRD_LON(g,l)           ) ! [OUT]
       enddo
       enddo

       ierr = fio_fclose(fid)
    else
       if( IO_L ) write(IO_FID_LOG,*) 'Invalid io_mode!'
       call PRC_abort
    endif

    return
  end subroutine MKIIMAP_input_src_hgrid

  !-----------------------------------------------------------------------------
  subroutine MKIIMAP_output_iimap
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_RGN_level, &
       PRC_RGN_local, &
       PRC_RGN_l2r
    use mod_fio, only: &
       FIO_output
    implicit none

    integer, allocatable :: prc_tab(:)

    character(len=H_MID) :: desc = 'IIMAP FILE'

    integer :: ierr
    !---------------------------------------------------------------------------

    ! reset fio parameters
    allocate( prc_tab(PRC_RGN_local) )
    prc_tab(1:PRC_RGN_local) = PRC_RGN_l2r(1:PRC_RGN_local)-1

    ierr = fio_put_commoninfo( IO_SPLIT_FILE,  & ! [IN]
                               IO_BIG_ENDIAN,  & ! [IN]
                               IO_ICOSAHEDRON, & ! [IN]
                               ADM_glevel,     & ! [IN]
                               PRC_RGN_level,  & ! [IN]
                               PRC_RGN_local,  & ! [IN]
                               prc_tab         ) ! [IN]

    deallocate(prc_tab)



    write(IO_FID_LOG,*) 'Output iimap. basename = ', trim(OUT_BASENAME)

    if ( OUT_io_mode == 'ADVANCED' ) then

       call FIO_output( iimap(:,:,:,I_rgnid), OUT_BASENAME, desc, '',     & ! [IN]
                       'iimap_rgnid', 'iimap(region id)', '',             & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( iimap(:,:,:,I_g1),    OUT_BASENAME, desc, '',     & ! [IN]
                       'iimap_g1',    'iimap(grid id 1)', '',             & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( iimap(:,:,:,I_g2),    OUT_BASENAME, desc, '',     & ! [IN]
                       'iimap_g2',    'iimap(grid id 2)', '',             & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( iimap(:,:,:,I_g3),    OUT_BASENAME, desc, '',     & ! [IN]
                       'iimap_g3',    'iimap(grid id 3)', '',             & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( iimap(:,:,:,I_w1),    OUT_BASENAME, desc, '',     & ! [IN]
                       'iimap_w1',    'iimap(grid weight 1)', '',         & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( iimap(:,:,:,I_w2),    OUT_BASENAME, desc, '',     & ! [IN]
                       'iimap_w2',    'iimap(grid weight 2)', '',         & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( iimap(:,:,:,I_w3),    OUT_BASENAME, desc, '',     & ! [IN]
                       'iimap_w3',    'iimap(grid weight 3)', '',         & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]

    else
       if( IO_L ) write(IO_FID_LOG,*) 'Invalid io_mode!'
       call PRC_abort
    endif

    return
  end subroutine MKIIMAP_output_iimap

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_mkiimap
