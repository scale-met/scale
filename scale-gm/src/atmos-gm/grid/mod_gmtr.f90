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

  use mod_adm, only: &
     ADM_TI,      &
     ADM_TJ,      &
     ADM_AI,      &
     ADM_AIJ,     &
     ADM_AJ,      &
     ADM_KNONE,   &
     ADM_lall,    &
     ADM_lall_pl, &
     ADM_gall,    &
     ADM_gall_pl
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
  integer,  public, parameter :: GMTR_p_nmax = 8

  integer,  public, parameter :: GMTR_p_AREA  = 1
  integer,  public, parameter :: GMTR_p_RAREA = 2
  integer,  public, parameter :: GMTR_p_IX    = 3
  integer,  public, parameter :: GMTR_p_IY    = 4
  integer,  public, parameter :: GMTR_p_IZ    = 5
  integer,  public, parameter :: GMTR_p_JX    = 6
  integer,  public, parameter :: GMTR_p_JY    = 7
  integer,  public, parameter :: GMTR_p_JZ    = 8

  integer,  public, parameter :: GMTR_t_nmax = 5

  integer,  public, parameter :: GMTR_t_AREA  = 1
  integer,  public, parameter :: GMTR_t_RAREA = 2
  integer,  public, parameter :: GMTR_t_W1    = 3
  integer,  public, parameter :: GMTR_t_W2    = 4
  integer,  public, parameter :: GMTR_t_W3    = 5

  integer,  public, parameter :: GMTR_a_nmax    = 12
  integer,  public, parameter :: GMTR_a_nmax_pl = 18

  integer,  public, parameter :: GMTR_a_HNX  = 1
  integer,  public, parameter :: GMTR_a_HNY  = 2
  integer,  public, parameter :: GMTR_a_HNZ  = 3
  integer,  public, parameter :: GMTR_a_HTX  = 4
  integer,  public, parameter :: GMTR_a_HTY  = 5
  integer,  public, parameter :: GMTR_a_HTZ  = 6
  integer,  public, parameter :: GMTR_a_TNX  = 7
  integer,  public, parameter :: GMTR_a_TNY  = 8
  integer,  public, parameter :: GMTR_a_TNZ  = 9
  integer,  public, parameter :: GMTR_a_TTX  = 10
  integer,  public, parameter :: GMTR_a_TTY  = 11
  integer,  public, parameter :: GMTR_a_TTZ  = 12

  integer,  public, parameter :: GMTR_a_TN2X = 13
  integer,  public, parameter :: GMTR_a_TN2Y = 14
  integer,  public, parameter :: GMTR_a_TN2Z = 15
  integer,  public, parameter :: GMTR_a_TT2X = 16
  integer,  public, parameter :: GMTR_a_TT2Y = 17
  integer,  public, parameter :: GMTR_a_TT2Z = 18

#ifdef _FIXEDINDEX_
  real(RP), public              :: GMTR_p   (ADM_gall   ,ADM_KNONE,ADM_lall   ,              GMTR_p_nmax   )
  real(RP), public              :: GMTR_p_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_p_nmax   )
  real(RP), public              :: GMTR_t   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,GMTR_t_nmax   )
  real(RP), public              :: GMTR_t_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_t_nmax   )
  real(RP), public              :: GMTR_a   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_AI:ADM_AJ,GMTR_a_nmax   )
  real(RP), public              :: GMTR_a_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_a_nmax_pl)

  real(RP), public              :: GMTR_area   (ADM_gall   ,ADM_lall   )
  real(RP), public              :: GMTR_area_pl(ADM_gall_pl,ADM_lall_pl)
#else
  real(RP), public, allocatable :: GMTR_p   (:,:,:,:)   ! geometrics for the cell point
  real(RP), public, allocatable :: GMTR_p_pl(:,:,:,:)
  real(RP), public, allocatable :: GMTR_t   (:,:,:,:,:) ! geometrics for the cell vertex
  real(RP), public, allocatable :: GMTR_t_pl(:,:,:,:)
  real(RP), public, allocatable :: GMTR_a   (:,:,:,:,:) ! geometrics for the cell arc
  real(RP), public, allocatable :: GMTR_a_pl(:,:,:,:)

  real(RP), public, allocatable :: GMTR_area   (:,:)    ! control area of the cell
  real(RP), public, allocatable :: GMTR_area_pl(:,:)
#endif

  character(len=H_SHORT), public :: GMTR_polygon_type = 'ON_SPHERE'
                                                      ! 'ON_SPHERE' triangle is fit to the sphere
                                                      ! 'ON_PLANE'  triangle is treated as 2D

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: GMTR_p_setup
  private :: GMTR_t_setup
  private :: GMTR_a_setup
  private :: GMTR_output_metrics

  private :: GMTR_TNvec

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: GMTR_fname   = ''
  character(len=H_SHORT), private :: GMTR_io_mode = 'ADVANCED'

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine GMTR_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_rscale, &
       GRD_x,      &
       GRD_x_pl,   &
       GRD_xt,     &
       GRD_xt_pl,  &
       GRD_s,      &
       GRD_s_pl
    implicit none

    namelist / GMTRPARAM / &
       GMTR_polygon_type, &
       GMTR_io_mode,      &
       GMTR_fname

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[gmtr]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=GMTRPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** GMTRPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist GMTRPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=GMTRPARAM)



#ifndef _FIXEDINDEX_
    allocate( GMTR_p   (ADM_gall   ,ADM_KNONE,ADM_lall   ,              GMTR_p_nmax   ) )
    allocate( GMTR_p_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_p_nmax   ) )
    allocate( GMTR_t   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,GMTR_t_nmax   ) )
    allocate( GMTR_t_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_t_nmax   ) )
    allocate( GMTR_a   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_AI:ADM_AJ,GMTR_a_nmax   ) )
    allocate( GMTR_a_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_a_nmax_pl) )

    allocate( GMTR_area    (ADM_gall,   ADM_lall   ) )
    allocate( GMTR_area_pl (ADM_gall_pl,ADM_lall_pl) )
#endif



    !--- calc geometrical information for cell point
    call GMTR_p_setup( GRD_x (:,:,:,:),   GRD_x_pl (:,:,:,:), & ! [IN]
                       GRD_xt(:,:,:,:,:), GRD_xt_pl(:,:,:,:), & ! [IN]
                       GRD_s (:,:,:,:),   GRD_s_pl (:,:,:,:), & ! [IN]
                       GMTR_p(:,:,:,:),   GMTR_p_pl(:,:,:,:), & ! [OUT]
                       GRD_rscale                             ) ! [IN]

    ! fill HALO
    call COMM_data_transfer( GMTR_p, GMTR_p_pl )

    ! for simple use
    GMTR_area   (:,:) = GMTR_p   (:,ADM_KNONE,:,GMTR_p_AREA)
    GMTR_area_pl(:,:) = GMTR_p_pl(:,ADM_KNONE,:,GMTR_p_AREA)

    !--- calc geometrical information for cell vertex (triangle)
    call GMTR_t_setup( GRD_x (:,:,:,:),   GRD_x_pl (:,:,:,:), & ! [IN]
                       GRD_xt(:,:,:,:,:), GRD_xt_pl(:,:,:,:), & ! [IN]
                       GMTR_t(:,:,:,:,:), GMTR_t_pl(:,:,:,:), & ! [OUT]
                       GRD_rscale                             ) ! [IN]

    !--- calc geometrical information for cell arc
    call GMTR_a_setup( GRD_x (:,:,:,:),   GRD_x_pl (:,:,:,:), & ! [IN]
                       GRD_xt(:,:,:,:,:), GRD_xt_pl(:,:,:,:), & ! [IN]
                       GMTR_a(:,:,:,:,:), GMTR_a_pl(:,:,:,:), & ! [OUT]
                       GRD_rscale                             ) ! [IN]

    if ( GMTR_fname /= "" ) then
       call GMTR_output_metrics( GMTR_fname )
    endif

    return
  end subroutine GMTR_setup

  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell point
  subroutine GMTR_p_setup( &
       GRD_x,  GRD_x_pl,  &
       GRD_xt, GRD_xt_pl, &
       GRD_s,  GRD_s_pl,  &
       GMTR_p, GMTR_p_pl, &
       GRD_rscale         )
    use mod_adm, only: &
       ADM_nxyz,     &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_vlink,    &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl
    use mod_grd, only: &
       GRD_XDIR,     &
       GRD_YDIR,     &
       GRD_ZDIR,     &
       GRD_LON,      &
       GRD_grid_type
    use scale_vector, only: &
       VECTR_triangle,      &
       VECTR_triangle_plane
    implicit none

    real(RP), intent(in)  :: GRD_x    (ADM_gall   ,ADM_KNONE,ADM_lall   ,              ADM_nxyz)
    real(RP), intent(in)  :: GRD_x_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
    real(RP), intent(in)  :: GRD_s    (ADM_gall   ,ADM_KNONE,ADM_lall   ,              2)
    real(RP), intent(in)  :: GRD_s_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              2)
    real(RP), intent(out) :: GMTR_p   (ADM_gall   ,ADM_KNONE,ADM_lall   ,GMTR_p_nmax)
    real(RP), intent(out) :: GMTR_p_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,GMTR_p_nmax)
    real(RP), intent(in)  :: GRD_rscale

    real(RP) :: wk   (ADM_nxyz,0:7,ADM_gall)
    real(RP) :: wk_pl(ADM_nxyz,0:ADM_vlink+1)

    real(RP) :: area
    real(RP) :: cos_lambda, sin_lambda

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: i, j, k0, l, d, v, n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** setup metrics for hexagonal/pentagonal mesh'

    k0 = ADM_KNONE

    GMTR_p   (:,:,:,:)   = 0.0_RP
    GMTR_p_pl(:,:,:,:)   = 0.0_RP

    do l = 1, ADM_lall
       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij     = suf(i  ,j  )
          ip1j   = suf(i+1,j  )
          ip1jp1 = suf(i+1,j+1)
          ijp1   = suf(i  ,j+1)
          im1j   = suf(i-1,j  )
          im1jm1 = suf(i-1,j-1)
          ijm1   = suf(i  ,j-1)

          !--- prepare 1 center and 6 vertices
          do d = 1, ADM_nxyz
             wk(d,0,ij) = GRD_x(ij,k0,l,d)

             wk(d,1,ij) = GRD_xt(ijm1  ,k0,l,ADM_TJ,d)
             wk(d,2,ij) = GRD_xt(ij    ,k0,l,ADM_TI,d)
             wk(d,3,ij) = GRD_xt(ij    ,k0,l,ADM_TJ,d)
             wk(d,4,ij) = GRD_xt(im1j  ,k0,l,ADM_TI,d)
             wk(d,5,ij) = GRD_xt(im1jm1,k0,l,ADM_TJ,d)
             wk(d,6,ij) = GRD_xt(im1jm1,k0,l,ADM_TI,d)
             wk(d,7,ij) = wk(d,1,ij)
          enddo
       enddo ! i loop
       enddo ! j loop

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,6,suf(ADM_gmin,ADM_gmin)) = wk(:,1,suf(ADM_gmin,ADM_gmin))
          wk(:,7,suf(ADM_gmin,ADM_gmin)) = wk(:,1,suf(ADM_gmin,ADM_gmin))
       endif

       !--- calc control area
       if ( GRD_grid_type == 'ON_PLANE' ) then
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = suf(i,j)

             area = 0.0_RP
             do v = 1, 6
                area = area + VECTR_triangle_plane( wk(:,0,ij), wk(:,v,ij), wk(:,v+1,ij) )
             enddo

             GMTR_p(ij,k0,l,GMTR_p_AREA)  = area
             GMTR_p(ij,k0,l,GMTR_p_RAREA) = 1.0_RP / GMTR_p(ij,k0,l,GMTR_p_AREA)

          enddo ! i loop
          enddo ! j loop
       else
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = suf(i,j)

             wk(:,:,ij) = wk(:,:,ij) / GRD_rscale

             area = 0.0_RP
             do v = 1, 6
                area = area + VECTR_triangle( wk(:,0,ij), wk(:,v,ij), wk(:,v+1,ij), GMTR_polygon_type, GRD_rscale )
             enddo

             GMTR_p(ij,k0,l,GMTR_p_AREA)  = area
             GMTR_p(ij,k0,l,GMTR_p_RAREA) = 1.0_RP / GMTR_p(ij,k0,l,GMTR_p_AREA)

          enddo ! i loop
          enddo ! j loop
       endif

       !--- calc coefficient between xyz <-> latlon
       if ( GRD_grid_type == 'ON_PLANE' ) then
          GMTR_p(:,k0,l,GMTR_p_IX) = 1.0_RP
          GMTR_p(:,k0,l,GMTR_p_IY) = 0.0_RP
          GMTR_p(:,k0,l,GMTR_p_IZ) = 0.0_RP
          GMTR_p(:,k0,l,GMTR_p_JX) = 0.0_RP
          GMTR_p(:,k0,l,GMTR_p_JY) = 1.0_RP
          GMTR_p(:,k0,l,GMTR_p_JZ) = 0.0_RP
       else
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = suf(i,j)

             sin_lambda = sin( GRD_s(ij,k0,l,GRD_LON) )
             cos_lambda = cos( GRD_s(ij,k0,l,GRD_LON) )

             GMTR_p(ij,k0,l,GMTR_p_IX) = -sin_lambda
             GMTR_p(ij,k0,l,GMTR_p_IY) =  cos_lambda
             GMTR_p(ij,k0,l,GMTR_p_IZ) = 0.0_RP
             GMTR_p(ij,k0,l,GMTR_p_JX) = -( GRD_x(ij,k0,l,GRD_ZDIR) * cos_lambda ) / GRD_rscale
             GMTR_p(ij,k0,l,GMTR_p_JY) = -( GRD_x(ij,k0,l,GRD_ZDIR) * sin_lambda ) / GRD_rscale
             GMTR_p(ij,k0,l,GMTR_p_JZ) =  ( GRD_x(ij,k0,l,GRD_XDIR) * cos_lambda &
                                          + GRD_x(ij,k0,l,GRD_YDIR) * sin_lambda ) / GRD_rscale
          enddo ! i loop
          enddo ! j loop
       endif
    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
          !--- prepare 1 center and * vertices
          do d = 1, ADM_nxyz
             wk_pl(d,0) = GRD_x_pl(n,k0,l,d)
             do v = 1, ADM_vlink ! (ICO=5)
                wk_pl(d,v) = GRD_xt_pl(v+1,k0,l,d)
             enddo
             wk_pl(d,ADM_vlink+1) = wk_pl(d,1)
          enddo

          wk_pl(:,:) = wk_pl(:,:) / GRD_rscale

          !--- calc control area
          area = 0.0_RP
          do v = 1, ADM_vlink ! (ICO=5)
             area = area + VECTR_triangle( wk_pl(:,0), wk_pl(:,v), wk_pl(:,v+1), GMTR_polygon_type, GRD_rscale )
          enddo

          GMTR_p_pl(n,k0,l,GMTR_p_AREA)  = area
          GMTR_p_pl(n,k0,l,GMTR_p_RAREA) = 1.0_RP / GMTR_p_pl(n,k0,l,GMTR_p_AREA)

          !--- calc coefficient between xyz <-> latlon
          sin_lambda = sin( GRD_s_pl(n,k0,l,GRD_LON) )
          cos_lambda = cos( GRD_s_pl(n,k0,l,GRD_LON) )

          GMTR_p_pl(n,k0,l,GMTR_p_IX) = -sin_lambda
          GMTR_p_pl(n,k0,l,GMTR_p_IY) =  cos_lambda
          GMTR_p_pl(n,k0,l,GMTR_p_IZ) = 0.0_RP
          GMTR_p_pl(n,k0,l,GMTR_p_JX) = -( GRD_x_pl(n,k0,l,GRD_ZDIR) * cos_lambda ) / GRD_rscale
          GMTR_p_pl(n,k0,l,GMTR_p_JY) = -( GRD_x_pl(n,k0,l,GRD_ZDIR) * sin_lambda ) / GRD_rscale
          GMTR_p_pl(n,k0,l,GMTR_p_JZ) =  ( GRD_x_pl(n,k0,l,GRD_XDIR) * cos_lambda &
                                         + GRD_x_pl(n,k0,l,GRD_YDIR) * sin_lambda ) / GRD_rscale
       enddo ! l loop
    endif

    return
  end subroutine GMTR_p_setup

  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell vertex (triangle)
  subroutine GMTR_t_setup( &
       GRD_x,  GRD_x_pl,  &
       GRD_xt, GRD_xt_pl, &
       GMTR_t, GMTR_t_pl, &
       GRD_rscale         )
    use mod_adm, only: &
       ADM_nxyz,     &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_grd, only: &
       GRD_grid_type
    use scale_vector, only: &
       VECTR_triangle,      &
       VECTR_triangle_plane
    implicit none

    real(RP), intent(in)  :: GRD_x    (ADM_gall   ,ADM_KNONE,ADM_lall   ,              ADM_nxyz)
    real(RP), intent(in)  :: GRD_x_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
    real(RP), intent(out) :: GMTR_t   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,GMTR_t_nmax)
    real(RP), intent(out) :: GMTR_t_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_t_nmax)
    real(RP), intent(in)  :: GRD_rscale

    real(RP) :: wk   (ADM_nxyz,0:3,ADM_gall,ADM_TI:ADM_TJ)
    real(RP) :: wk_pl(ADM_nxyz,0:3)

    real(RP) :: area, area1, area2, area3

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1

    integer  :: i, j, k0, l, d, v, n, t
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** setup metrics for triangle mesh'

    k0 = ADM_KNONE

    GMTR_t   (:,:,:,:,:) = 0.0_RP
    GMTR_t_pl(:,:,:,:)   = 0.0_RP

    do l = 1,ADM_lall
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij     = suf(i  ,j  )
          ip1j   = suf(i+1,j  )
          ip1jp1 = suf(i+1,j+1)
          ijp1   = suf(i  ,j+1)

          !--- prepare 1 center and 3 vertices for 2 triangles
          do d = 1, ADM_nxyz
             wk(d,0,ij,ADM_TI) = GRD_xt(ij,k0,l,ADM_TI,d)

             wk(d,1,ij,ADM_TI) = GRD_x(ij    ,k0,l,d)
             wk(d,2,ij,ADM_TI) = GRD_x(ip1j  ,k0,l,d)
             wk(d,3,ij,ADM_TI) = GRD_x(ip1jp1,k0,l,d)

             wk(d,0,ij,ADM_TJ) = GRD_xt(ij,k0,l,ADM_TJ,d)

             wk(d,1,ij,ADM_TJ) = GRD_x(ij    ,k0,l,d)
             wk(d,2,ij,ADM_TJ) = GRD_x(ip1jp1,k0,l,d)
             wk(d,3,ij,ADM_TJ) = GRD_x(ijp1  ,k0,l,d)
          enddo
       enddo
       enddo

       !--- treat unused triangle
       wk(:,:,suf(ADM_gmax,ADM_gmin-1),ADM_TI) = wk(:,:,suf(ADM_gmax,ADM_gmin-1),ADM_TJ)
       wk(:,:,suf(ADM_gmin-1,ADM_gmax),ADM_TJ) = wk(:,:,suf(ADM_gmin-1,ADM_gmax),ADM_TI)

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,:,suf(ADM_gmin-1,ADM_gmin-1),ADM_TI) = wk(:,:,suf(ADM_gmin,ADM_gmin-1),ADM_TJ)
       endif

       if ( GRD_grid_type == 'ON_PLANE' ) then
          do t = ADM_TI,ADM_TJ
          do j = ADM_gmin-1, ADM_gmax
          do i = ADM_gmin-1, ADM_gmax
             ij = suf(i,j)

             area1 = VECTR_triangle_plane( wk(:,0,ij,t), wk(:,2,ij,t), wk(:,3,ij,t) )
             area2 = VECTR_triangle_plane( wk(:,0,ij,t), wk(:,3,ij,t), wk(:,1,ij,t) )
             area3 = VECTR_triangle_plane( wk(:,0,ij,t), wk(:,1,ij,t), wk(:,2,ij,t) )

             area = area1 + area2 + area3

             GMTR_t(ij,k0,l,t,GMTR_t_AREA)  = area
             GMTR_t(ij,k0,l,t,GMTR_t_RAREA) = 1.0_RP / area

             GMTR_t(ij,k0,l,t,GMTR_t_W1)    = area1 / area
             GMTR_t(ij,k0,l,t,GMTR_t_W2)    = area2 / area
             GMTR_t(ij,k0,l,t,GMTR_t_W3)    = area3 / area
          enddo
          enddo
          enddo
       else
          do t = ADM_TI,ADM_TJ
          do j = ADM_gmin-1, ADM_gmax
          do i = ADM_gmin-1, ADM_gmax
             ij = suf(i,j)

             wk(:,:,ij,t) = wk(:,:,ij,t) / GRD_rscale

             area1 = VECTR_triangle( wk(:,0,ij,t), wk(:,2,ij,t), wk(:,3,ij,t), GMTR_polygon_type, GRD_rscale )
             area2 = VECTR_triangle( wk(:,0,ij,t), wk(:,3,ij,t), wk(:,1,ij,t), GMTR_polygon_type, GRD_rscale )
             area3 = VECTR_triangle( wk(:,0,ij,t), wk(:,1,ij,t), wk(:,2,ij,t), GMTR_polygon_type, GRD_rscale )

             area = area1 + area2 + area3

             GMTR_t(ij,k0,l,t,GMTR_t_AREA)  = area
             GMTR_t(ij,k0,l,t,GMTR_t_RAREA) = 1.0_RP / area

             GMTR_t(ij,k0,l,t,GMTR_t_W1)    = area1 / area
             GMTR_t(ij,k0,l,t,GMTR_t_W2)    = area2 / area
             GMTR_t(ij,k0,l,t,GMTR_t_W3)    = area3 / area
          enddo
          enddo
          enddo
       endif

    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1,ADM_lall_pl
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             do d = 1, ADM_nxyz
                wk_pl(d,0) = GRD_xt_pl(ij,k0,l,d)

                wk_pl(d,1) = GRD_x_pl(n   ,k0,l,d)
                wk_pl(d,2) = GRD_x_pl(ij  ,k0,l,d)
                wk_pl(d,3) = GRD_x_pl(ijp1,k0,l,d)
             enddo

             wk_pl(:,:) = wk_pl(:,:) / GRD_rscale

             area1 = VECTR_triangle( wk_pl(:,0), wk_pl(:,2), wk_pl(:,3), GMTR_polygon_type, GRD_rscale )
             area2 = VECTR_triangle( wk_pl(:,0), wk_pl(:,3), wk_pl(:,1), GMTR_polygon_type, GRD_rscale )
             area3 = VECTR_triangle( wk_pl(:,0), wk_pl(:,1), wk_pl(:,2), GMTR_polygon_type, GRD_rscale )

             area = area1 + area2 + area3

             GMTR_t_pl(ij,k0,l,GMTR_t_AREA)  = area
             GMTR_t_pl(ij,k0,l,GMTR_t_RAREA) = 1.0_RP / area

             GMTR_t_pl(ij,k0,l,GMTR_t_W1)    = area1 / area
             GMTR_t_pl(ij,k0,l,GMTR_t_W2)    = area2 / area
             GMTR_t_pl(ij,k0,l,GMTR_t_W3)    = area3 / area
          enddo
       enddo
    endif

    return
  end subroutine GMTR_t_setup

  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell arc
  subroutine GMTR_a_setup( &
       GRD_x,  GRD_x_pl,  &
       GRD_xt, GRD_xt_pl, &
       GMTR_a, GMTR_a_pl, &
       GRD_rscale         )
    use mod_adm, only: &
       ADM_nxyz,     &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_grd, only: &
       GRD_grid_type
    implicit none

    real(RP), intent(in)  :: GRD_x    (ADM_gall   ,ADM_KNONE,ADM_lall   ,              ADM_nxyz)
    real(RP), intent(in)  :: GRD_x_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
    real(RP), intent(out) :: GMTR_a   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_AI:ADM_AJ,GMTR_a_nmax   )
    real(RP), intent(out) :: GMTR_a_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_a_nmax_pl)
    real(RP), intent(in)  :: GRD_rscale

    real(RP) :: wk   (ADM_nxyz,2,ADM_gall)
    real(RP) :: wk_pl(ADM_nxyz,2)

    real(RP) :: Tvec(3), Nvec(3)

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1

    integer  :: i, j, k0, l, d, v, n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** setup metrics for cell arcs'

    k0 = ADM_KNONE

    GMTR_a   (:,:,:,:,:) = 0.0_RP
    GMTR_a_pl(:,:,:,:)   = 0.0_RP

    !--- Triangle
    do l = 1, ADM_lall

       !--- AI
       do j = ADM_gmin-1, ADM_gmax+1
       do i = ADM_gmin-1, ADM_gmax
          ij   = suf(i  ,j  )
          ip1j = suf(i+1,j  )

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_x(ij  ,k0,l,d)
             wk(d,2,ij) = GRD_x(ip1j,k0,l,d)
          enddo
       enddo
       enddo

       ! treat arc of unused triangle
       wk(:,1,suf(ADM_gmax  ,ADM_gmin-1)) = GRD_x(suf(ADM_gmax  ,ADM_gmin-1),k0,l,:)
       wk(:,2,suf(ADM_gmax  ,ADM_gmin-1)) = GRD_x(suf(ADM_gmax  ,ADM_gmin  ),k0,l,:)
       wk(:,1,suf(ADM_gmin-1,ADM_gmax+1)) = GRD_x(suf(ADM_gmin  ,ADM_gmax+1),k0,l,:)
       wk(:,2,suf(ADM_gmin-1,ADM_gmax+1)) = GRD_x(suf(ADM_gmin  ,ADM_gmax  ),k0,l,:)

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,1,suf(ADM_gmin-1,ADM_gmin-1)) = GRD_x(suf(ADM_gmin  ,ADM_gmin-1),k0,l,:)
          wk(:,2,suf(ADM_gmin-1,ADM_gmin-1)) = GRD_x(suf(ADM_gmin+1,ADM_gmin  ),k0,l,:)
       endif

       do j = ADM_gmin-1, ADM_gmax+1
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_TNX) = Nvec(1)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_TNY) = Nvec(2)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_TNZ) = Nvec(3)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_TTX) = Tvec(1)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_TTY) = Tvec(2)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_TTZ) = Tvec(3)
       enddo
       enddo

       !--- AIJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij     = suf(i  ,j  )
          ip1jp1 = suf(i+1,j+1)

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_x(ij    ,k0,l,d)
             wk(d,2,ij) = GRD_x(ip1jp1,k0,l,d)
          enddo
       enddo
       enddo

       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_TNX) = Nvec(1)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_TNY) = Nvec(2)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_TNZ) = Nvec(3)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_TTX) = Tvec(1)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_TTY) = Tvec(2)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_TTZ) = Tvec(3)
       enddo
       enddo

       !--- AJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax+1
          ij   = suf(i  ,j  )
          ijp1 = suf(i  ,j+1)

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_x(ij  ,k0,l,d)
             wk(d,2,ij) = GRD_x(ijp1,k0,l,d)
          enddo
       enddo
       enddo

       ! treat arc of unused triangle
       wk(:,1,suf(ADM_gmax+1,ADM_gmin-1)) = GRD_x(suf(ADM_gmax+1,ADM_gmin),k0,l,:)
       wk(:,2,suf(ADM_gmax+1,ADM_gmin-1)) = GRD_x(suf(ADM_gmax  ,ADM_gmin),k0,l,:)
       wk(:,1,suf(ADM_gmin-1,ADM_gmax  )) = GRD_x(suf(ADM_gmin-1,ADM_gmax),k0,l,:)
       wk(:,2,suf(ADM_gmin-1,ADM_gmax  )) = GRD_x(suf(ADM_gmin  ,ADM_gmax),k0,l,:)

       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax+1
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_TNX) = Nvec(1)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_TNY) = Nvec(2)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_TNZ) = Nvec(3)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_TTX) = Tvec(1)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_TTY) = Tvec(2)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_TTZ) = Tvec(3)
       enddo
       enddo

    enddo ! l loop

    !--- Hexagon
    do l = 1, ADM_lall

       !--- AI
       do j = ADM_gmin,   ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij   = suf(i  ,j  )
          ijm1 = suf(i  ,j-1)

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_xt(ij  ,k0,l,ADM_TI,d)
             wk(d,2,ij) = GRD_xt(ijm1,k0,l,ADM_TJ,d)
          enddo
       enddo
       enddo

       do j = ADM_gmin,   ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_HNX) = Nvec(1)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_HNY) = Nvec(2)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_HNZ) = Nvec(3)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_HTX) = Tvec(1)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_HTY) = Tvec(2)
          GMTR_a(ij,k0,l,ADM_AI,GMTR_a_HTZ) = Tvec(3)
       enddo
       enddo

       !--- AIJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij   = suf(i  ,j  )

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_xt(ij  ,k0,l,ADM_TJ,d)
             wk(d,2,ij) = GRD_xt(ij  ,k0,l,ADM_TI,d)
          enddo
       enddo
       enddo

       ! treat arc of unused hexagon
       wk(:,1,suf(ADM_gmax  ,ADM_gmin-1)) = GRD_xt(suf(ADM_gmax  ,ADM_gmin-1),k0,l,ADM_TJ,:)
       wk(:,2,suf(ADM_gmax  ,ADM_gmin-1)) = GRD_xt(suf(ADM_gmax  ,ADM_gmin  ),k0,l,ADM_TI,:)
       wk(:,1,suf(ADM_gmin-1,ADM_gmax  )) = GRD_xt(suf(ADM_gmin  ,ADM_gmax  ),k0,l,ADM_TJ,:)
       wk(:,2,suf(ADM_gmin-1,ADM_gmax  )) = GRD_xt(suf(ADM_gmin-1,ADM_gmax  ),k0,l,ADM_TI,:)

       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_HNX) = Nvec(1)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_HNY) = Nvec(2)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_HNZ) = Nvec(3)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_HTX) = Tvec(1)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_HTY) = Tvec(2)
          GMTR_a(ij,k0,l,ADM_AIJ,GMTR_a_HTZ) = Tvec(3)
       enddo
       enddo

       !--- AJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin,   ADM_gmax
          ij   = suf(i  ,j  )
          im1j = suf(i-1,j  )

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_xt(im1j,k0,l,ADM_TI,d)
             wk(d,2,ij) = GRD_xt(ij  ,k0,l,ADM_TJ,d)
          enddo
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,1,suf(ADM_gmin  ,ADM_gmin-1)) = GRD_xt(suf(ADM_gmin  ,ADM_gmin  ),k0,l,ADM_TI,:)
          wk(:,2,suf(ADM_gmin  ,ADM_gmin-1)) = GRD_xt(suf(ADM_gmin  ,ADM_gmin-1),k0,l,ADM_TJ,:)
       endif

       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin,   ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_HNX) = Nvec(1)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_HNY) = Nvec(2)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_HNZ) = Nvec(3)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_HTX) = Tvec(1)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_HTY) = Tvec(2)
          GMTR_a(ij,k0,l,ADM_AJ,GMTR_a_HTZ) = Tvec(3)
       enddo
       enddo

    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl

          !--- Triangle (arc 1)
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij = v

             do d = 1, ADM_nxyz
                wk_pl(d,1) = GRD_x_pl(n ,k0,l,d)
                wk_pl(d,2) = GRD_x_pl(ij,k0,l,d)
             enddo

             call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                              wk_pl(:,1), wk_pl(:,2),                      & ! [IN]
                              GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

             GMTR_a_pl(ij,k0,l,GMTR_a_TNX) = Nvec(1)
             GMTR_a_pl(ij,k0,l,GMTR_a_TNY) = Nvec(2)
             GMTR_a_pl(ij,k0,l,GMTR_a_TNZ) = Nvec(3)
             GMTR_a_pl(ij,k0,l,GMTR_a_TTX) = Tvec(1)
             GMTR_a_pl(ij,k0,l,GMTR_a_TTY) = Tvec(2)
             GMTR_a_pl(ij,k0,l,GMTR_a_TTZ) = Tvec(3)
          enddo

          !--- Triangle (arc 2)
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v+1
             if ( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             do d = 1, ADM_nxyz
                wk_pl(d,1) = GRD_x_pl(ij  ,k0,l,d)
                wk_pl(d,2) = GRD_x_pl(ijp1,k0,l,d)
             enddo

             call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                              wk_pl(:,1), wk_pl(:,2),                      & ! [IN]
                              GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

             GMTR_a_pl(ij,k0,l,GMTR_a_TN2X) = Nvec(1)
             GMTR_a_pl(ij,k0,l,GMTR_a_TN2Y) = Nvec(2)
             GMTR_a_pl(ij,k0,l,GMTR_a_TN2Z) = Nvec(3)
             GMTR_a_pl(ij,k0,l,GMTR_a_TT2X) = Tvec(1)
             GMTR_a_pl(ij,k0,l,GMTR_a_TT2Y) = Tvec(2)
             GMTR_a_pl(ij,k0,l,GMTR_a_TT2Z) = Tvec(3)
          enddo

          !--- Hexagon
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v-1
             if ( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl

             do d = 1, ADM_nxyz
                wk_pl(d,1) = GRD_xt_pl(ijm1,k0,l,d)
                wk_pl(d,2) = GRD_xt_pl(ij  ,k0,l,d)
             enddo

             call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                              wk_pl(:,1), wk_pl(:,2),                      & ! [IN]
                              GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

             GMTR_a_pl(ij,k0,l,GMTR_a_HNX) = Nvec(1)
             GMTR_a_pl(ij,k0,l,GMTR_a_HNY) = Nvec(2)
             GMTR_a_pl(ij,k0,l,GMTR_a_HNZ) = Nvec(3)
             GMTR_a_pl(ij,k0,l,GMTR_a_HTX) = Tvec(1)
             GMTR_a_pl(ij,k0,l,GMTR_a_HTY) = Tvec(2)
             GMTR_a_pl(ij,k0,l,GMTR_a_HTZ) = Tvec(3)
          enddo

       enddo
    endif

    return
  end subroutine GMTR_a_setup

  !-----------------------------------------------------------------------------
  subroutine GMTR_TNvec( &
       vT,           &
       vN,           &
       vFrom,        &
       vTo,          &
       grid_type,    &
       polygon_type, &
       radius        )
    use scale_vector, only: &
       VECTR_dot,   &
       VECTR_cross, &
       VECTR_abs,   &
       VECTR_angle
    implicit none

    real(RP),         intent(out) :: vT   (3)     ! tangential vector
    real(RP),         intent(out) :: vN   (3)     ! normal     vector
    real(RP),         intent(in)  :: vFrom(3)
    real(RP),         intent(in)  :: vTo  (3)
    character(len=*), intent(in)  :: grid_type    ! ON_SPHERE or ON_PLANE
    character(len=*), intent(in)  :: polygon_type ! ON_SPHERE or ON_PLANE
    real(RP),         intent(in)  :: radius

    real(RP), parameter :: o(3) = 0.0_RP

    real(RP) :: angle, length
    real(RP) :: distance
    !---------------------------------------------------------------------------

    ! calculate tangential vector
    vT(:) = vTo(:) - vFrom(:)

    if ( grid_type == 'ON_PLANE' ) then ! treat as point on the plane

       ! calculate normal vector
       vN(1) = -vT(2)
       vN(2) =  vT(1)
       vN(3) = 0.0_RP

    elseif( grid_type == 'ON_SPHERE' ) then ! treat as point on the sphere

       if ( polygon_type == 'ON_PLANE' ) then ! length of a line

          call VECTR_dot( distance, vFrom(:), vTo(:), vFrom(:), vTo(:) )
          distance = sqrt( distance )

       elseif( polygon_type == 'ON_SPHERE' ) then ! length of a geodesic line ( angle * radius )

          call VECTR_angle( angle, vFrom(:), o(:), vTo(:) )
          distance = angle * radius

       endif

       call VECTR_abs( length, vT(:) )
       vT(:) = vT(:) * distance / length

       ! calculate normal vector
       call VECTR_cross( vN(:), o(:), vFrom(:), o(:), vTo(:) )

       call VECTR_abs( length, vN(:) )
       vN(:) = vN(:) * distance / length

    endif

    return
  end subroutine GMTR_TNvec

  !-----------------------------------------------------------------------------
  subroutine GMTR_output_metrics( &
       basename )
    use mod_io_param, only: &
       IO_REAL8, &
       IO_REAL4
    use scale_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_have_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_fio, only: &
       FIO_output
    implicit none

    character(len=*), intent(in) :: basename

    character(len=H_MID)  :: desc = 'Metrics info'

    real(RP) :: tmp    (ADM_gall   ,ADM_KNONE,ADM_lall   ,2)
    real(RP) :: tmp_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,2)

    real(RP) :: tmp2   (ADM_gall   ,54,ADM_lall   ,1)
    real(RP) :: tmp2_pl(ADM_gall_pl,54,ADM_lall_pl,1)

    integer :: rgnid
    integer,  parameter :: I_rgn  = 1
    integer,  parameter :: I_grid = 2

    integer :: dtype
    integer :: g, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    if    ( RP == SP ) then
       dtype = IO_REAL4
    elseif( RP == DP ) then
       dtype = IO_REAL8
    endif

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)
       do g = 1, ADM_gall
          tmp(g,k0,l,I_rgn ) = real(rgnid,kind=RP)
          tmp(g,k0,l,I_grid) = real(g    ,kind=RP)
       enddo
    enddo

    if ( ADM_have_pl ) Then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl
          tmp_pl(g,k0,l,I_rgn ) = real(-l,kind=RP)
          tmp_pl(g,k0,l,I_grid) = real(g ,kind=RP)
       enddo
       enddo
    endif

    call COMM_data_transfer( tmp, tmp_pl )

    do l = 1, ADM_lall
    do g = 1, ADM_gall
          tmp2(g, 1: 8,l,1) = abs( GMTR_p(g,k0,l        ,:) )
          tmp2(g, 9:13,l,1) = abs( GMTR_t(g,k0,l,ADM_TI ,:) )
          tmp2(g,14:18,l,1) = abs( GMTR_t(g,k0,l,ADM_TJ ,:) )
          tmp2(g,19:30,l,1) = abs( GMTR_a(g,k0,l,ADM_AI ,:) )
          tmp2(g,31:42,l,1) = abs( GMTR_a(g,k0,l,ADM_AIJ,:) )
          tmp2(g,43:54,l,1) = abs( GMTR_a(g,k0,l,ADM_AJ ,:) )
    enddo
    enddo

    if ( ADM_have_pl ) Then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl
          tmp2_pl(g, 1: 8,l,1) = abs( GMTR_p_pl(g,k0,l,:) )
          tmp2_pl(g, 9:13,l,1) = abs( GMTR_t_pl(g,k0,l,:) )
          tmp2_pl(g,14:18,l,1) = abs( GMTR_t_pl(g,k0,l,:) )
          tmp2_pl(g,19:30,l,1) = abs( GMTR_a_pl(g,k0,l,1:GMTR_a_nmax) )
          tmp2_pl(g,31:42,l,1) = abs( GMTR_a_pl(g,k0,l,1:GMTR_a_nmax) )
          tmp2_pl(g,43:54,l,1) = abs( GMTR_a_pl(g,k0,l,1:GMTR_a_nmax) )
       enddo
       enddo
    endif

    call COMM_data_transfer( tmp2, tmp2_pl )

    if ( GMTR_io_mode == 'ADVANCED' ) then
       call FIO_output( tmp(:,:,:,I_rgn),  basename, desc, "",          & ! [IN]
                        "rgn", "region number", "",                     & ! [IN]
                        "NIL", dtype, "ZSSFC1", 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( tmp(:,:,:,I_grid), basename, desc, "",          & ! [IN]
                        "grid", "grid number", "",                      & ! [IN]
                        "NIL", dtype, "ZSSFC1", 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( tmp2(:,:,:,1),     basename, desc, "",          & ! [IN]
                        "gmtrmetrics", "gmtr metrics", "",              & ! [IN]
                        "", dtype, "LAYERNM", 1, 54, 1, 0.0_DP, 0.0_DP  ) ! [IN]
    endif

    return
  end subroutine GMTR_output_metrics

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
