!-------------------------------------------------------------------------------
!> Module geometrics
!!
!! @par Description
!!         In this module, the geometrics of the icosahedral grid such as area are calculated
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_gmtr
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_LOG_FID
  use mod_adm, only: &
     ADM_TI,      &
     ADM_TJ,      &
     ADM_AI,      &
     ADM_AIJ,     &
     ADM_AJ,      &
     ADM_lall,    &
     ADM_lall_pl, &
     ADM_gall,    &
     ADM_gall_pl, &
     ADM_kall,    &
     ADM_KNONE
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GMTR_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: GMTR_P_nmax_var = 10

  integer, public, parameter :: GMTR_P_AREA  = 1
  integer, public, parameter :: GMTR_P_RAREA = 2
  integer, public, parameter :: GMTR_P_IX    = 3
  integer, public, parameter :: GMTR_P_IY    = 4
  integer, public, parameter :: GMTR_P_IZ    = 5
  integer, public, parameter :: GMTR_P_JX    = 6
  integer, public, parameter :: GMTR_P_JY    = 7
  integer, public, parameter :: GMTR_P_JZ    = 8
  integer, public, parameter :: GMTR_P_LAT   = 9
  integer, public, parameter :: GMTR_P_LON   = 10

  integer, public, parameter :: GMTR_T_nmax_var = 7

  integer, public, parameter :: GMTR_T_AREA  = 1
  integer, public, parameter :: GMTR_T_RAREA = 2
  integer, public, parameter :: GMTR_T_W1    = 3
  integer, public, parameter :: GMTR_T_W2    = 4
  integer, public, parameter :: GMTR_T_W3    = 5
  integer, public, parameter :: GMTR_T_LAT   = 6
  integer, public, parameter :: GMTR_T_LON   = 7

  integer, public, parameter :: GMTR_A_nmax_var    = 12
  integer, public, parameter :: GMTR_A_nmax_var_pl = 18

  integer, public, parameter :: GMTR_A_HNX  = 1
  integer, public, parameter :: GMTR_A_HNY  = 2
  integer, public, parameter :: GMTR_A_HNZ  = 3
  integer, public, parameter :: GMTR_A_HTX  = 4
  integer, public, parameter :: GMTR_A_HTY  = 5
  integer, public, parameter :: GMTR_A_HTZ  = 6
  integer, public, parameter :: GMTR_A_TNX  = 7
  integer, public, parameter :: GMTR_A_TNY  = 8
  integer, public, parameter :: GMTR_A_TNZ  = 9
  integer, public, parameter :: GMTR_A_TTX  = 10
  integer, public, parameter :: GMTR_A_TTY  = 11
  integer, public, parameter :: GMTR_A_TTZ  = 12

  integer, public, parameter :: GMTR_A_TN2X = 13
  integer, public, parameter :: GMTR_A_TN2Y = 14
  integer, public, parameter :: GMTR_A_TN2Z = 15
  integer, public, parameter :: GMTR_A_TT2X = 16
  integer, public, parameter :: GMTR_A_TT2Y = 17
  integer, public, parameter :: GMTR_A_TT2Z = 18

#ifdef _FIXEDINDEX_
  real(RP), public              :: GMTR_P_var   (ADM_gall   ,ADM_KNONE,ADM_lall   ,              GMTR_P_nmax_var   )
  real(RP), public              :: GMTR_P_var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_P_nmax_var   )
  real(RP), public              :: GMTR_T_var   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,GMTR_T_nmax_var   )
  real(RP), public              :: GMTR_T_var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_T_nmax_var   )
  real(RP), public              :: GMTR_A_var   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_AI:ADM_AJ,GMTR_A_nmax_var   )
  real(RP), public              :: GMTR_A_var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_A_nmax_var_pl)

  real(RP), public              :: GMTR_area    (ADM_gall   ,ADM_lall   )
  real(RP), public              :: GMTR_area_pl (ADM_gall_pl,ADM_lall_pl)
  real(RP), public              :: GMTR_lat     (ADM_gall   ,ADM_lall   )
  real(RP), public              :: GMTR_lat_pl  (ADM_gall_pl,ADM_lall_pl)
  real(RP), public              :: GMTR_lon     (ADM_gall   ,ADM_lall   )
  real(RP), public              :: GMTR_lon_pl  (ADM_gall_pl,ADM_lall_pl)
#else
  real(RP), public, allocatable :: GMTR_P_var   (:,:,:,:)   ! geometrics for the cell point
  real(RP), public, allocatable :: GMTR_P_var_pl(:,:,:,:)
  real(RP), public, allocatable :: GMTR_T_var   (:,:,:,:,:) ! geometrics for the cell vertex
  real(RP), public, allocatable :: GMTR_T_var_pl(:,:,:,:)
  real(RP), public, allocatable :: GMTR_A_var   (:,:,:,:,:) ! geometrics for the cell arc
  real(RP), public, allocatable :: GMTR_A_var_pl(:,:,:,:)

  real(RP), public, allocatable :: GMTR_area    (:,:)       ! control area of the cell
  real(RP), public, allocatable :: GMTR_area_pl (:,:)
  real(RP), public, allocatable :: GMTR_lat     (:,:)       ! latitude  of the cell point
  real(RP), public, allocatable :: GMTR_lat_pl  (:,:)
  real(RP), public, allocatable :: GMTR_lon     (:,:)       ! longitude of the cell point
  real(RP), public, allocatable :: GMTR_lon_pl  (:,:)
#endif

  character(len=H_SHORT), public :: GMTR_polygon_type = 'ON_SPHERE'
                                                       ! 'ON_SPHERE' triangle is fit to the sphere
                                                       ! 'ON_PLANE'  triangle is treated as 2D

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: GMTR_calc_P
  private :: GMTR_calc_T
  private :: GMTR_calc_A
  private :: GMTR_output_metrics

  private :: mk_gmtrvec_on_plane
  private :: triangle_area_on_plane

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private :: GMTR_fname   = ''
  character(len=H_SHORT),     private :: GMTR_io_mode = 'LEGACY'

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine GMTR_setup
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_gmin,      &
       ADM_gmax
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    character(len=H_SHORT) :: polygon_type

    namelist / GMTRPARAM / &
       polygon_type

    integer :: ierr
    integer :: K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    polygon_type = GMTR_polygon_type

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[gmtr]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=GMTRPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** GMTRPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist GMTRPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist GMTRPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=GMTRPARAM)

    GMTR_polygon_type = polygon_type



#ifndef _FIXEDINDEX_
    allocate( GMTR_P_var   (ADM_gall   ,K0,ADM_lall   ,              GMTR_P_nmax_var   ) )
    allocate( GMTR_P_var_pl(ADM_gall_pl,K0,ADM_lall_pl,              GMTR_P_nmax_var   ) )
    allocate( GMTR_T_var   (ADM_gall   ,K0,ADM_lall   ,ADM_TI:ADM_TJ,GMTR_T_nmax_var   ) )
    allocate( GMTR_T_var_pl(ADM_gall_pl,K0,ADM_lall_pl,              GMTR_T_nmax_var   ) )
    allocate( GMTR_A_var   (ADM_gall   ,K0,ADM_lall   ,ADM_AI:ADM_AJ,GMTR_A_nmax_var   ) )
    allocate( GMTR_A_var_pl(ADM_gall_pl,K0,ADM_lall_pl,              GMTR_A_nmax_var_pl) )

    allocate( GMTR_area    (ADM_gall,   ADM_lall   ) )
    allocate( GMTR_area_pl (ADM_gall_pl,ADM_lall_pl) )
    allocate( GMTR_lat     (ADM_gall,   ADM_lall   ) )
    allocate( GMTR_lat_pl  (ADM_gall_pl,ADM_lall_pl) )
    allocate( GMTR_lon     (ADM_gall,   ADM_lall   ) )
    allocate( GMTR_lon_pl  (ADM_gall_pl,ADM_lall_pl) )
#endif
    GMTR_P_var   (:,:,:,:)   = 0.0_RP
    GMTR_P_var_pl(:,:,:,:)   = 0.0_RP
    GMTR_T_var   (:,:,:,:,:) = 0.0_RP
    GMTR_T_var_pl(:,:,:,:)   = 0.0_RP
    GMTR_A_var   (:,:,:,:,:) = 0.0_RP
    GMTR_A_var_pl(:,:,:,:)   = 0.0_RP



    !--- calc geometrical information for cell point
    call GMTR_calc_P

    !--- calc geometrical information for cell vertex (triangle)
    call GMTR_calc_T

    !--- calc geometrical information for cell arc
    call GMTR_calc_A



    !--- fill HALO
    call COMM_data_transfer( GMTR_P_var, GMTR_P_var_pl )

    GMTR_P_var(suf(ADM_gmax+1,ADM_gmin-1),:,:,:) = GMTR_P_var(suf(ADM_gmax+1,ADM_gmin),:,:,:)
    GMTR_P_var(suf(ADM_gmin-1,ADM_gmax+1),:,:,:) = GMTR_P_var(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    !--- for simple use
    GMTR_area   (:,:) = GMTR_P_var   (:,K0,:,GMTR_P_AREA)
    GMTR_area_pl(:,:) = GMTR_P_var_pl(:,K0,:,GMTR_P_AREA)
    GMTR_lat    (:,:) = GMTR_P_var   (:,K0,:,GMTR_P_LAT )
    GMTR_lat_pl (:,:) = GMTR_P_var_pl(:,K0,:,GMTR_P_LAT )
    GMTR_lon    (:,:) = GMTR_P_var   (:,K0,:,GMTR_P_LON )
    GMTR_lon_pl (:,:) = GMTR_P_var_pl(:,K0,:,GMTR_P_LON )

    if ( GMTR_fname /= "" ) then
       call GMTR_output_metrics( GMTR_fname )
    endif

    return
  end subroutine GMTR_setup

  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell point
  subroutine GMTR_calc_P
    use scale_vector, only: &
       VECTR_xyz2latlon, &
       VECTR_triangle
    use mod_adm, only: &
       ADM_prc_me,      &
       ADM_have_pl,     &
       ADM_prc_tab,     &
       ADM_rgn_vnum,    &
       ADM_W,           &
       ADM_nxyz,        &
       ADM_gmin,        &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GIoJm,       &
       ADM_GImJo,       &
       ADM_GImJm,       &
       ADM_gslf_pl,     &
       ADM_gmin_pl,     &
       ADM_vlink_nmax
    use mod_grd, only: &
       GRD_XDIR,      &
       GRD_YDIR,      &
       GRD_ZDIR,      &
       GRD_x,         &
       GRD_x_pl,      &
       GRD_xt,        &
       GRD_xt_pl,     &
       GRD_grid_type, &
       GRD_rscale
    implicit none

    real(RP) :: v   (ADM_nxyz,0:7,ADM_gall)
    real(RP) :: v_pl(ADM_nxyz,0:ADM_vlink_nmax+1)

    real(RP) :: area
    real(RP) :: cos_lam, sin_lam

    integer :: l, n, m
    integer :: rgnid, ij, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do n = 1, ADM_IooJoo_nmax
          ij = ADM_IooJoo(n,ADM_GIoJo)

          v(GRD_XDIR,0,ij) = GRD_x(ij,K0,l,GRD_XDIR)
          v(GRD_XDIR,1,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,3,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,4,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJo),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,5,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,6,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,7,ij) = v(GRD_XDIR,1,ij)

          v(GRD_YDIR,0,ij) = GRD_x(ij,K0,l,GRD_YDIR)
          v(GRD_YDIR,1,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,3,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,4,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJo),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,5,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,6,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,7,ij) = v(GRD_YDIR,1,ij)

          v(GRD_ZDIR,0,ij) = GRD_x(ij,K0,l,GRD_ZDIR)
          v(GRD_ZDIR,1,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,3,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,4,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJo),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,5,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,6,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,7,ij) = v(GRD_ZDIR,1,ij)
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          v(:,6,suf(ADM_gmin,ADM_gmin)) = v(:,1,suf(ADM_gmin,ADM_gmin))
          v(:,7,suf(ADM_gmin,ADM_gmin)) = v(:,1,suf(ADM_gmin,ADM_gmin))
       endif

       do n = 1, ADM_IooJoo_nmax
          ij = ADM_IooJoo(n,ADM_GIoJo)

          area = 0.0_RP
          if ( GRD_grid_type == 'ON_PLANE' ) then
             do m = 1, 6
                area = area + triangle_area_on_plane( v(:,0,ij), v(:,m,ij), v(:,m+1,ij) )
             enddo
          else
             do m = 1, 6
                area = area + VECTR_triangle( v(:,0,ij), v(:,m,ij), v(:,m+1,ij), GMTR_polygon_type, GRD_rscale )
             enddo
          endif

          GMTR_P_var(ij,K0,l,GMTR_P_AREA)  = area
          GMTR_P_var(ij,K0,l,GMTR_P_RAREA) = 1.0_RP / GMTR_P_var(ij,K0,l,GMTR_P_AREA)

          call VECTR_xyz2latlon( GRD_x     (ij,K0,l,GRD_XDIR),   & ! [IN]
                                 GRD_x     (ij,K0,l,GRD_YDIR),   & ! [IN]
                                 GRD_x     (ij,K0,l,GRD_ZDIR),   & ! [IN]
                                 GMTR_P_var(ij,K0,l,GMTR_P_LAT), & ! [OUT]
                                 GMTR_P_var(ij,K0,l,GMTR_P_LON)  ) ! [OUT]

          if ( GRD_grid_type == 'ON_PLANE' ) then

             GMTR_P_var(ij,K0,l,GMTR_P_IX) = 1.0_RP
             GMTR_P_var(ij,K0,l,GMTR_P_IY) = 0.0_RP
             GMTR_P_var(ij,K0,l,GMTR_P_IZ) = 0.0_RP
             GMTR_P_var(ij,K0,l,GMTR_P_JX) = 0.0_RP
             GMTR_P_var(ij,K0,l,GMTR_P_JY) = 1.0_RP
             GMTR_P_var(ij,K0,l,GMTR_P_JZ) = 0.0_RP

          else

             sin_lam = sin( GMTR_P_var(ij,K0,l,GMTR_P_LON) )
             cos_lam = cos( GMTR_P_var(ij,K0,l,GMTR_P_LON) )

             GMTR_P_var(ij,K0,l,GMTR_P_IX) = -sin_lam
             GMTR_P_var(ij,K0,l,GMTR_P_IY) =  cos_lam
             GMTR_P_var(ij,K0,l,GMTR_P_IZ) = 0.0_RP
             GMTR_P_var(ij,K0,l,GMTR_P_JX) = -( GRD_x(ij,K0,l,GRD_ZDIR) * cos_lam ) / GRD_rscale
             GMTR_P_var(ij,K0,l,GMTR_P_JY) = -( GRD_x(ij,K0,l,GRD_ZDIR) * sin_lam ) / GRD_rscale
             GMTR_P_var(ij,K0,l,GMTR_P_JZ) =  ( GRD_x(ij,K0,l,GRD_XDIR) * cos_lam &
                                              + GRD_x(ij,K0,l,GRD_YDIR) * sin_lam ) / GRD_rscale
          endif

       enddo ! ij loop
    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_GSLF_PL

       do l = 1,ADM_lall_pl

          v_pl(GRD_XDIR,0) = GRD_x_pl(n,K0,l,GRD_XDIR)
          v_pl(GRD_YDIR,0) = GRD_x_pl(n,K0,l,GRD_YDIR)
          v_pl(GRD_ZDIR,0) = GRD_x_pl(n,K0,l,GRD_ZDIR)
          do m = 1, ADM_vlink_nmax ! (ICO=5)
             v_pl(GRD_XDIR,m) = GRD_xt_pl(m+ADM_GMIN_PL-1,K0,l,GRD_XDIR)
             v_pl(GRD_YDIR,m) = GRD_xt_pl(m+ADM_GMIN_PL-1,K0,l,GRD_YDIR)
             v_pl(GRD_ZDIR,m) = GRD_xt_pl(m+ADM_GMIN_PL-1,K0,l,GRD_ZDIR)
          enddo
          v_pl(:,ADM_vlink_nmax+1) = v_pl(:,1)

          area = 0.0_RP
          do m = 1, ADM_vlink_nmax ! (ICO=5)
             area = area + VECTR_triangle( v_pl(:,0), v_pl(:,m), v_pl(:,m+1), GMTR_polygon_type, GRD_rscale )
          enddo

          GMTR_P_var_pl(n,K0,l,GMTR_P_AREA)  = area
          GMTR_P_var_pl(n,K0,l,GMTR_P_RAREA) = 1.0_RP / GMTR_P_var_pl(n,K0,l,GMTR_P_AREA)

          call VECTR_xyz2latlon( GRD_x_pl     (n,K0,l,GRD_XDIR),   & ! [IN]
                                 GRD_x_pl     (n,K0,l,GRD_YDIR),   & ! [IN]
                                 GRD_x_pl     (n,K0,l,GRD_ZDIR),   & ! [IN]
                                 GMTR_P_var_pl(n,K0,l,GMTR_P_LAT), & ! [OUT]
                                 GMTR_P_var_pl(n,K0,l,GMTR_P_LON)  ) ! [OUT]

          sin_lam = sin( GMTR_P_var_pl(n,K0,l,GMTR_P_LON) )
          cos_lam = cos( GMTR_P_var_pl(n,K0,l,GMTR_P_LON) )

          GMTR_P_var_pl(n,K0,l,GMTR_P_IX) = -sin_lam
          GMTR_P_var_pl(n,K0,l,GMTR_P_IY) =  cos_lam
          GMTR_P_var_pl(n,K0,l,GMTR_P_IZ) = 0.0_RP
          GMTR_P_var_pl(n,K0,l,GMTR_P_JX) = -( GRD_x_pl(n,K0,l,GRD_ZDIR) * cos_lam ) / GRD_rscale
          GMTR_P_var_pl(n,K0,l,GMTR_P_JY) = -( GRD_x_pl(n,K0,l,GRD_ZDIR) * sin_lam ) / GRD_rscale
          GMTR_P_var_pl(n,K0,l,GMTR_P_JZ) =  ( GRD_x_pl(n,K0,l,GRD_XDIR) * cos_lam &
                                             + GRD_x_pl(n,K0,l,GRD_YDIR) * sin_lam ) / GRD_rscale
       enddo ! l loop
    endif

    return
  end subroutine GMTR_calc_P

  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell vertex (triangle)
  subroutine GMTR_calc_T
    use scale_vector, only: &
       VECTR_xyz2latlon, &
       VECTR_triangle
    use mod_adm, only: &
       ADM_prc_me,      &
       ADM_have_pl,     &
       ADM_prc_tab,     &
       ADM_rgn_vnum,    &
       ADM_W,           &
       ADM_nxyz,        &
       ADM_gmin,        &
       ADM_gmax,        &
       ADM_ImoJmo_nmax, &
       ADM_ImoJmo,      &
       ADM_GIoJo,       &
       ADM_GIpJo,       &
       ADM_GIpJp,       &
       ADM_GIoJp,       &
       ADM_gslf_pl,     &
       ADM_gmin_pl,     &
       ADM_gmax_pl
    use mod_grd, only: &
       GRD_XDIR,      &
       GRD_YDIR,      &
       GRD_ZDIR,      &
       GRD_x,         &
       GRD_x_pl,      &
       GRD_xt,        &
       GRD_xt_pl,     &
       GRD_grid_type, &
       GRD_rscale
    implicit none

    real(RP) :: v   (ADM_nxyz,0:3,ADM_gall   ,ADM_TI:ADM_TJ)
    real(RP) :: v_pl(ADM_nxyz,0:3,ADM_gall_pl)

    real(RP) :: area, area1, area2, area3
    integer :: l, d, t, n
    integer :: rgnid, ij, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    do l = 1,ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          v(GRD_XDIR,0,ij,ADM_TI) = GRD_xt(ij,K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,1,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,3,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_XDIR)

          v(GRD_XDIR,0,ij,ADM_TJ) = GRD_xt(ij,K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,1,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_XDIR)
          v(GRD_XDIR,3,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJp),K0,l,GRD_XDIR)

          v(GRD_YDIR,0,ij,ADM_TI) = GRD_xt(ij,K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,1,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,3,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_YDIR)

          v(GRD_YDIR,0,ij,ADM_TJ) = GRD_xt(ij,K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,1,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_YDIR)
          v(GRD_YDIR,3,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJp),K0,l,GRD_YDIR)

          v(GRD_ZDIR,0,ij,ADM_TI) = GRD_xt(ij,K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,1,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,3,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_ZDIR)

          v(GRD_ZDIR,0,ij,ADM_TJ) = GRD_xt(ij,K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,1,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,3,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJp),K0,l,GRD_ZDIR)
       enddo

       !--- treat unused point
       v(:,:,suf(ADM_gmax,ADM_gmin-1),ADM_TI) = v(:,:,suf(ADM_gmax,ADM_gmin-1),ADM_TJ)
       v(:,:,suf(ADM_gmin-1,ADM_gmax),ADM_TJ) = v(:,:,suf(ADM_gmin-1,ADM_gmax),ADM_TI)

       !--- exception for the west
       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          v(:,:,suf(ADM_gmin-1,ADM_gmin-1),ADM_TI) = v(:,:,suf(ADM_gmin,ADM_gmin-1),ADM_TJ)
       endif

       do t = ADM_TI,ADM_TJ
       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          if ( GRD_grid_type == 'ON_PLANE' ) then
             area1 = triangle_area_on_plane( v(:,0,ij,t), v(:,2,ij,t), v(:,3,ij,t) )
             area2 = triangle_area_on_plane( v(:,0,ij,t), v(:,3,ij,t), v(:,1,ij,t) )
             area3 = triangle_area_on_plane( v(:,0,ij,t), v(:,1,ij,t), v(:,2,ij,t) )
          else
             area1 = VECTR_triangle( v(:,0,ij,t), v(:,2,ij,t), v(:,3,ij,t), GMTR_polygon_type, GRD_rscale )
             area2 = VECTR_triangle( v(:,0,ij,t), v(:,3,ij,t), v(:,1,ij,t), GMTR_polygon_type, GRD_rscale )
             area3 = VECTR_triangle( v(:,0,ij,t), v(:,1,ij,t), v(:,2,ij,t), GMTR_polygon_type, GRD_rscale )
          endif

          area = area1 + area2 + area3

          GMTR_T_var(ij,K0,l,t,GMTR_T_AREA)  = area
          GMTR_T_var(ij,K0,l,t,GMTR_T_RAREA) = 1.0_RP / area

          GMTR_T_var(ij,K0,l,t,GMTR_T_W1)    = area1 / area
          GMTR_T_var(ij,K0,l,t,GMTR_T_W2)    = area2 / area
          GMTR_T_var(ij,K0,l,t,GMTR_T_W3)    = area3 / area

          call VECTR_xyz2latlon( GRD_xt    (ij,K0,l,t,GRD_XDIR),   & ! [IN]
                                 GRD_xt    (ij,K0,l,t,GRD_YDIR),   & ! [IN]
                                 GRD_xt    (ij,K0,l,t,GRD_ZDIR),   & ! [IN]
                                 GMTR_T_var(ij,K0,l,t,GMTR_T_LAT), & ! [OUT]
                                 GMTR_T_var(ij,K0,l,t,GMTR_T_LON)  ) ! [OUT]
       enddo
       enddo

    enddo

    if ( ADM_have_pl ) then
       do l = 1,ADM_lall_pl

          do n = ADM_GMIN_PL, ADM_GMAX_PL-1
             do d = GRD_XDIR, GRD_ZDIR
                v_pl(d,0,n) = GRD_xt_pl(n,K0,l,d)
                v_pl(d,1,n) = GRD_x_pl(ADM_GSLF_PL,K0,l,d)
                v_pl(d,2,n) = GRD_x_pl(n  ,        K0,l,d)
                v_pl(d,3,n) = GRD_x_pl(n+1,        K0,l,d)
             enddo
          enddo

          n = ADM_GMAX_PL
          do d = GRD_XDIR, GRD_ZDIR
             v_pl(d,0,n) = GRD_xt_pl(n,K0,l,d)
             v_pl(d,1,n) = GRD_x_pl(ADM_GSLF_PL,K0,l,d)
             v_pl(d,2,n) = GRD_x_pl(n,          K0,l,d)
             v_pl(d,3,n) = GRD_x_pl(ADM_GMIN_PL,K0,l,d)
          enddo

          do n = ADM_GMIN_PL, ADM_GMAX_PL
             area1 = VECTR_triangle( v_pl(:,0,n), v_pl(:,2,n), v_pl(:,3,n), GMTR_polygon_type, GRD_rscale )
             area2 = VECTR_triangle( v_pl(:,0,n), v_pl(:,3,n), v_pl(:,1,n), GMTR_polygon_type, GRD_rscale )
             area3 = VECTR_triangle( v_pl(:,0,n), v_pl(:,1,n), v_pl(:,2,n), GMTR_polygon_type, GRD_rscale )

             area = area1 + area2 + area3

             GMTR_T_var_pl(n,K0,l,GMTR_T_AREA)  = area
             GMTR_T_var_pl(n,K0,l,GMTR_T_RAREA) = 1.0_RP / area

             GMTR_T_var_pl(n,K0,l,GMTR_T_W1)    = area1 / area
             GMTR_T_var_pl(n,K0,l,GMTR_T_W2)    = area2 / area
             GMTR_T_var_pl(n,K0,l,GMTR_T_W3)    = area3 / area

             call VECTR_xyz2latlon( GRD_xt_pl    (n,K0,l,GRD_XDIR),   & ! [IN]
                                    GRD_xt_pl    (n,K0,l,GRD_YDIR),   & ! [IN]
                                    GRD_xt_pl    (n,K0,l,GRD_ZDIR),   & ! [IN]
                                    GMTR_T_var_pl(n,K0,l,GMTR_T_LAT), & ! [OUT]
                                    GMTR_T_var_pl(n,K0,l,GMTR_T_LON)  ) ! [OUT]
          enddo

       enddo
    endif

    return
  end subroutine GMTR_calc_T

  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell arc
  subroutine GMTR_calc_A
    use scale_vector, only: &
       GMTR_calc_vector
    use mod_adm, only: &
       ADM_prc_me,      &
       ADM_have_pl,     &
       ADM_prc_tab,     &
       ADM_rgn_vnum,    &
       ADM_W,           &
       ADM_nxyz,        &
       ADM_gmin,        &
       ADM_gmax,        &
       ADM_ImpJmo_nmax, &
       ADM_ImoJoo_nmax, &
       ADM_IooJmo_nmax, &
       ADM_ImoJmp_nmax, &
       ADM_ImoJmo_nmax, &
       ADM_ImpJmo,      &
       ADM_ImoJoo,      &
       ADM_IooJmo,      &
       ADM_ImoJmo,      &
       ADM_ImoJmp,      &
       ADM_GIoJo,       &
       ADM_GIpJo,       &
       ADM_GIpJp,       &
       ADM_GIoJp,       &
       ADM_GIoJm,       &
       ADM_GImJo,       &
       ADM_gslf_pl,     &
       ADM_gmin_pl,     &
       ADM_gmax_pl
    use mod_grd, only: &
       GRD_XDIR,      &
       GRD_YDIR,      &
       GRD_ZDIR,      &
       GRD_x,         &
       GRD_x_pl,      &
       GRD_xt,        &
       GRD_xt_pl,     &
       GRD_grid_type, &
       GRD_rscale
    implicit none

    real(RP) :: v   (ADM_nxyz,2,ADM_gall   )
    real(RP) :: v_pl(ADM_nxyz,2,ADM_gall_pl)

    real(RP) :: tvec(ADM_nxyz)
    real(RP) :: nvec(ADM_nxyz)

    integer :: ij, K0, l, d
    integer :: rgnid
    integer :: n
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    !--- Triangle
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       !--- AI
       do n = 1, ADM_ImoJmp_nmax
          ij = ADM_ImoJmp(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIpJo),K0,l,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIpJo),K0,l,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIpJo),K0,l,GRD_ZDIR)
       enddo

       do d = GRD_XDIR, GRD_ZDIR
          !--- execetion for the south.
          v(d,1,suf(ADM_gmax,ADM_gmin-1)) = GRD_x(suf(ADM_gmax,ADM_gmin-1),K0,l,d)
          v(d,2,suf(ADM_gmax,ADM_gmin-1)) = GRD_x(suf(ADM_gmax,ADM_gmin),  K0,l,d)

          !--- execetion for the south.
          v(d,1,suf(ADM_gmin-1,ADM_gmax+1)) = GRD_x(suf(ADM_gmin,ADM_gmax+1),K0,l,d)
          v(d,2,suf(ADM_gmin-1,ADM_gmax+1)) = GRD_x(suf(ADM_gmin,ADM_gmax),  K0,l,d)

          !--- exception for the west
          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
             v(d,1,suf(ADM_gmin-1,ADM_gmin-1)) = GRD_x(suf(ADM_gmin,ADM_gmin-1),K0,l,d)
             v(d,2,suf(ADM_gmin-1,ADM_gmin-1)) = GRD_x(suf(ADM_gmin+1,ADM_gmin),K0,l,d)
          endif
       enddo

       do n = 1, ADM_ImoJmp_nmax
          ij = ADM_ImoJmp(n,ADM_GIoJo)

          if ( GRD_grid_type == 'ON_PLANE' ) then
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call GMTR_calc_vector( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TNZ) = nvec(3)
       enddo

       !--- AIJ
       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_ZDIR)
       enddo

       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          if ( GRD_grid_type == 'ON_PLANE' ) then
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call GMTR_calc_vector( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TNZ) = nvec(3)
       enddo

       !--- AJ
       do n = 1, ADM_ImpJmo_nmax
          ij = ADM_ImpJmo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJp),K0,l,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJp),K0,l,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJp),K0,l,GRD_ZDIR)
       enddo

       do d = GRD_XDIR, GRD_ZDIR
          !--- execetion for the south.
          v(d,1,suf(ADM_gmax+1,ADM_gmin-1)) = GRD_x(suf(ADM_gmax+1,ADM_gmin),K0,l,d)
          v(d,2,suf(ADM_gmax+1,ADM_gmin-1)) = GRD_x(suf(ADM_gmax,ADM_gmin),  K0,l,d)

          !--- execetion for the north.
          v(d,1,suf(ADM_gmin-1,ADM_gmax)) = GRD_x(suf(ADM_gmin-1,ADM_gmax),K0,l,d)
          v(d,2,suf(ADM_gmin-1,ADM_gmax)) = GRD_x(suf(ADM_gmin,ADM_gmax),  K0,l,d)
       enddo

       do n = 1, ADM_ImpJmo_nmax
          ij = ADM_ImpJmo(n,ADM_GIoJo)

          if ( GRD_grid_type == 'ON_PLANE' ) then
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call GMTR_calc_vector( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TNZ) = nvec(3)
       enddo
    enddo

    !
    ! --- Hexagon
    !
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       !--- AI
       do n = 1, ADM_ImoJoo_nmax
          ij = ADM_ImoJoo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_ZDIR)
       enddo

       do n = 1, ADM_ImoJoo_nmax
          ij = ADM_ImoJoo(n,ADM_GIoJo)

          if ( GRD_grid_type == 'ON_PLANE' ) then
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call GMTR_calc_vector( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HNZ) = nvec(3)
       enddo

       !--- AIJ
       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_ZDIR)
       enddo

       do d = GRD_XDIR, GRD_ZDIR
          !--- execetion for the south.
          v(d,1,suf(ADM_gmax,ADM_gmin-1)) = GRD_xt(suf(ADM_gmax,ADM_gmin-1),K0,l,ADM_TJ,d)
          v(d,2,suf(ADM_gmax,ADM_gmin-1)) = GRD_xt(suf(ADM_gmax,ADM_gmin),  K0,l,ADM_TI,d)

          !--- execetion for the north.
          v(d,1,suf(ADM_gmin-1,ADM_gmax)) = GRD_xt(suf(ADM_gmin,ADM_gmax),  K0,l,ADM_TJ,d)
          v(d,2,suf(ADM_gmin-1,ADM_gmax)) = GRD_xt(suf(ADM_gmin-1,ADM_gmax),K0,l,ADM_TI,d)
       enddo

       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          if ( GRD_grid_type == 'ON_PLANE' ) then
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call GMTR_calc_vector( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HNZ) = nvec(3)
       enddo

       !--- AJ
       do n = 1, ADM_IooJmo_nmax
          ij = ADM_IooJmo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_xt(ADM_IooJmo(n,ADM_GImJo),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_xt(ADM_IooJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_xt(ADM_IooJmo(n,ADM_GImJo),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_xt(ADM_IooJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_xt(ADM_IooJmo(n,ADM_GImJo),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_xt(ADM_IooJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_ZDIR)
       enddo

       !--- exception for the west
       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          do d = GRD_XDIR, GRD_ZDIR
             v(d,1,suf(ADM_gmin,ADM_gmin-1)) = GRD_xt(suf(ADM_gmin,ADM_gmin),  K0,l,ADM_TI,d)
             v(d,2,suf(ADM_gmin,ADM_gmin-1)) = GRD_xt(suf(ADM_gmin,ADM_gmin-1),K0,l,ADM_TJ,d)
          enddo
       endif

       do n = 1, ADM_IooJmo_nmax
          ij = ADM_IooJmo(n,ADM_GIoJo)

          if ( GRD_grid_type == 'ON_PLANE' ) then
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call GMTR_calc_vector( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HNZ) = nvec(3)
       enddo

    enddo ! l loop

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl

          do ij = ADM_GMIN_PL, ADM_GMAX_PL
             do d = GRD_XDIR, GRD_ZDIR
                v_pl(d,1,ij) = GRD_x_pl(ADM_GSLF_PL,K0,l,d)
                v_pl(d,2,ij) = GRD_x_pl(ij         ,K0,l,d)
             enddo

             call GMTR_calc_vector( v_pl(:,1,ij), v_pl(:,2,ij), tvec(:), nvec(:), &
                                   GMTR_polygon_type, GRD_rscale                 )

             GMTR_A_var_pl(ij,K0,l,GMTR_A_TTX:GMTR_A_TTZ) = tvec(1:3)
             GMTR_A_var_pl(ij,K0,l,GMTR_A_TNX:GMTR_A_TNZ) = nvec(1:3)
          enddo

          do ij = ADM_GMIN_PL, ADM_GMAX_PL-1
             do d = GRD_XDIR, GRD_ZDIR
                v_pl(d,1,ij) = GRD_x_pl(ij,  K0,l,d)
                v_pl(d,2,ij) = GRD_x_pl(ij+1,K0,l,d)
             enddo
          enddo
          do d = GRD_XDIR, GRD_ZDIR
             v_pl(d,1,ADM_GMAX_PL) = GRD_x_pl(ADM_GMAX_PL,K0,l,d)
             v_pl(d,2,ADM_GMAX_PL) = GRD_x_pl(ADM_GMIN_PL,K0,l,d)
          enddo

          do ij = ADM_GMIN_PL, ADM_GMAX_PL
             call GMTR_calc_vector( v_pl(:,1,ij), v_pl(:,2,ij), tvec(:), nvec(:), &
                                   GMTR_polygon_type, GRD_rscale                 )

             GMTR_A_var_pl(ij,K0,l,GMTR_A_TT2X:GMTR_A_TT2Z) = tvec(1:3)
             GMTR_A_var_pl(ij,K0,l,GMTR_A_TN2X:GMTR_A_TN2Z) = nvec(1:3)
          enddo

          do d = GRD_XDIR, GRD_ZDIR
             v_pl(d,1,ADM_GMIN_PL) = GRD_xt_pl(ADM_GMAX_PL,K0,l,d)
             v_pl(d,2,ADM_GMIN_PL) = GRD_xt_pl(ADM_GMIN_PL,K0,l,d)
          enddo
          do ij = ADM_GMIN_PL+1, ADM_GMAX_PL
             do d = GRD_XDIR, GRD_ZDIR
                v_pl(d,1,ij) = GRD_xt_pl(ij-1,K0,l,d)
                v_pl(d,2,ij) = GRD_xt_pl(ij,  K0,l,d)
             enddo
          enddo

          do ij = ADM_GMIN_PL, ADM_GMAX_PL
             call GMTR_calc_vector( v_pl(:,1,ij), v_pl(:,2,ij), tvec(:), nvec(:), &
                                   GMTR_polygon_type, GRD_rscale                 )

             GMTR_A_var_pl(ij,K0,l,GMTR_A_HTX:GMTR_A_HTZ) = tvec(1:3)
             GMTR_A_var_pl(ij,K0,l,GMTR_A_HNX:GMTR_A_HNZ) = nvec(1:3)
          enddo

       enddo
    endif

    return
  end subroutine GMTR_calc_A

  !-----------------------------------------------------------------------------
  subroutine GMTR_output_metrics( &
       basename )
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_prc_tab,   &
       ADM_have_pl, &
       ADM_prc_me
    use mod_fio, only: &
       FIO_output, &
       FIO_REAL8
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    character(len=*), intent(in) :: basename

    character(len=H_LONG) :: fname
    character(len=H_MID)  :: desc = 'Metrics info'

    real(RP) :: tmp   (ADM_gall   ,ADM_KNONE,ADM_lall   ,2)
    real(RP) :: tmp_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,2)

    integer :: rgnid
    integer, parameter :: I_rgn  = 1
    integer, parameter :: I_grid = 2

    integer :: fid
    integer :: g, l, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)
       do g = 1, ADM_gall
          tmp(g,K0,l,I_rgn ) = real(rgnid,kind=RP)
          tmp(g,K0,l,I_grid) = real(g    ,kind=RP)
       enddo
    enddo

    if ( ADM_have_pl ) Then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl
          tmp_pl(g,K0,l,I_rgn ) = real(-l,kind=RP)
          tmp_pl(g,K0,l,I_grid) = real(g ,kind=RP)
       enddo
       enddo
    endif

    call COMM_data_transfer( tmp, tmp_pl )

    if ( GMTR_io_mode == 'ADVANCED' ) then

       call FIO_output( GMTR_P_var(:,:,:,GMTR_P_AREA),                     &
                        basename, desc, "",                                &
                        "area", "control area", "",                        &
                        "m2", FIO_REAL8, "ZSSFC1", 1, 1, 1, 0.0_RP, 0.0_RP     )
       call FIO_output( GMTR_P_var(:,:,:,GMTR_P_LAT),                      &
                        basename, desc, "",                                &
                        "lat", "latitude", "",                             &
                        "radian", FIO_REAL8, "ZSSFC1", 1, 1, 1, 0.0_RP, 0.0_RP )
       call FIO_output( GMTR_P_var(:,:,:,GMTR_P_LON),                      &
                        basename, desc, "",                                &
                        "lon", "longitude", "",                            &
                        "radian", FIO_REAL8, "ZSSFC1", 1, 1, 1, 0.0_RP, 0.0_RP )
       call FIO_output( tmp(:,:,:,I_rgn),                                  &
                        basename, desc, "",                                &
                        "rgn", "region number", "",                        &
                        "NIL", FIO_REAL8, "ZSSFC1", 1, 1, 1, 0.0_RP, 0.0_RP    )
       call FIO_output( tmp(:,:,:,I_grid),                                 &
                        basename, desc, "",                                &
                        "grid", "grid number", "",                         &
                        "NIL", FIO_REAL8, "ZSSFC1", 1, 1, 1, 0.0_RP, 0.0_RP    )

    elseif( GMTR_io_mode == 'LEGACY' ) then

       do l = 1, ADM_lall
          rgnid = ADM_prc_tab(l,ADM_prc_me)
          call IO_make_idstr(fname,trim(basename),'rgn',rgnid-1)

          fid = IO_get_available_fid()
          open( unit   = fid,           &
                file   = trim(fname),   &
                form   = 'unformatted', &
                access = 'direct',      &
                recl   = ADM_gall*8     )

             write(fid,rec=1) GMTR_P_var(:,K0,l,GMTR_P_AREA)
             write(fid,rec=2) GMTR_P_var(:,K0,l,GMTR_P_LAT )
             write(fid,rec=3) GMTR_P_var(:,K0,l,GMTR_P_LON )
             write(fid,rec=4) tmp       (:,K0,l,I_rgn      )
             write(fid,rec=5) tmp       (:,K0,l,I_grid     )

          close(fid)
       enddo

    else
       write(ADM_LOG_FID,*) 'Invalid io_mode!'
       call ADM_proc_stop
    endif

    return
  end subroutine GMTR_output_metrics

  !-----------------------------------------------------------------------------
  !> calc vector on plane
  subroutine mk_gmtrvec_on_plane( vFrom, vTo, vT, vN )
    implicit none

    real(RP), intent(in)  :: vFrom(3), vTo(3)
    real(RP), intent(out) :: vT(3),    vN(3)
    !---------------------------------------------------------------------------

    vT(:) = vTo(:) - vFrom(:)

    vN(1) = -vT(2)
    vN(2) =  vT(1)
    vN(3) =   0.0_RP

    return
  end subroutine mk_gmtrvec_on_plane

  !-----------------------------------------------------------------------------
  !> calc triangle area on plane
  !> @return area
  function triangle_area_on_plane( a, b, c ) result(area)
    implicit none

    real(RP), intent(in) :: a(3), b(3), c(3)
    real(RP)             :: area
    !
    real(RP) :: a2b(3), a2c(3)
    real(RP) :: len_a2b, len_a2c
    real(RP) :: prd
    !---------------------------------------------------------------------------

    a2b(:)  = b(:) - a(:)
    a2c(:)  = c(:) - a(:)
    len_a2b = a2b(1)*a2b(1) + a2b(3)*a2b(3) + a2b(3)*a2b(3) ! |a->b|**2
    len_a2c = a2c(1)*a2c(1) + a2c(3)*a2c(3) + a2c(3)*a2c(3) ! |a->c|**2
    prd     = a2b(1)*a2c(1) + a2b(2)*a2c(2) + a2b(3)*a2c(3) ! (a->b)*(a->c)

    area    = 0.5_RP * sqrt( len_a2b * len_a2c - prd*prd )

  end function triangle_area_on_plane

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_gmtr
