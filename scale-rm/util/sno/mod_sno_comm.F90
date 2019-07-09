!-------------------------------------------------------------------------------
!> Module SNO (RM)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          SCALE NetCDF Operator (SNO)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_sno_comm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SNO_comm_globalaxis
  public :: SNO_comm_globalvars

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  ! 1D
  private :: gather_x
  private :: gather_y
  ! 2D
  private :: gather_xy
  private :: gather_nx
  private :: gather_ny
  ! 3D
  private :: gather_nxy
  private :: gather_xyn

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNO_comm_globalaxis( &
       ismaster,      &
       output_single, &
       nprocs_x_out,  &
       nprocs_y_out,  &
       hinfo,         &
       naxis,         &
       ainfo,         &
       ainfo_all,     &
       bcast          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo
    implicit none

    logical,          intent(in)    :: ismaster                              ! master rank process?
    logical,          intent(in)    :: output_single                         ! output single file when using MPI?
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology
    type(commoninfo), intent(in)    :: hinfo                                 ! common information
    integer,          intent(in)    :: naxis                                 ! number of axis variables
    type(axisinfo),   intent(in)    :: ainfo(:)                              ! axis information (input)
    type(axisinfo),   intent(out)   :: ainfo_all(:)                          ! axis information (output)

    logical,          intent(in), optional :: bcast                          ! use broadcast?

    integer  :: D1, D2, D3
    integer  :: GD1, GD2, GD3
    integer  :: xhalo, yhalo
    integer  :: n
    integer  :: ierr

    logical  :: bcast_
    !---------------------------------------------------------------------------

    if ( output_single ) then

       xhalo = hinfo%halosize(2)
       yhalo = hinfo%halosize(3)

       if ( present( bcast ) ) then
          bcast_ = .true.
       else
          bcast_ = .false.
       endif

       do n = 1, naxis
          ! set axis information
          ainfo_all(n)%varname     = ainfo(n)%varname
          ainfo_all(n)%description = ainfo(n)%description
          ainfo_all(n)%units       = ainfo(n)%units
          ainfo_all(n)%datatype    = ainfo(n)%datatype
          ainfo_all(n)%dim_rank    = ainfo(n)%dim_rank
          ainfo_all(n)%dim_name(:) = ainfo(n)%dim_name(:)
          ainfo_all(n)%transpose   = ainfo(n)%transpose
          ainfo_all(n)%regrid      = ainfo(n)%regrid
          ainfo_all(n)%has_bounds  = ainfo(n)%has_bounds
          ainfo_all(n)%is_bounds   = ainfo(n)%is_bounds

          select case( ainfo_all(n)%dim_rank )
          case( 1 )

             D1 = size( ainfo(n)%AXIS_1d(:), 1 )

             select case( ainfo_all(n)%varname )
             case('x')

                if ( ismaster ) then
                   GD1 = D1 * nprocs_x_out
                endif
                call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                call gather_x( ismaster,               &
                               bcast_,                 &
                               ainfo_all(n)%datatype,  &
                               nprocs_x_out,           &
                               nprocs_y_out,           &
                               D1, 0, 0, 0,            &
                               ainfo(n)%AXIS_1d(:),    &
                               ainfo_all(n)%AXIS_1d(:) )

             case('xh','lon')

                if ( ismaster ) then
                   GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                endif
                call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                call gather_x( ismaster,               &
                               bcast_,                 &
                               ainfo_all(n)%datatype,  &
                               nprocs_x_out,           &
                               nprocs_y_out,           &
                               D1, 1, 0, 0,            &
                               ainfo(n)%AXIS_1d(:),    &
                               ainfo_all(n)%AXIS_1d(:) )

             case('y')

                if ( ismaster ) then
                   GD1 = D1 * nprocs_y_out
                endif
                call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                call gather_y( ismaster,               &
                               bcast_,                 &
                               ainfo_all(n)%datatype,  &
                               nprocs_x_out,           &
                               nprocs_y_out,           &
                               D1, 0, 0, 0,            &
                               ainfo(n)%AXIS_1d(:),    &
                               ainfo_all(n)%AXIS_1d(:) )

             case('yh','lat')

                if ( ismaster ) then
                   GD1 = ( D1 - 1 ) * nprocs_y_out + 1
                endif
                call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                call gather_y( ismaster,               &
                               bcast_,                 &
                               ainfo_all(n)%datatype,  &
                               nprocs_x_out,           &
                               nprocs_y_out,           &
                               D1, 1, 0, 0,            &
                               ainfo(n)%AXIS_1d(:),    &
                               ainfo_all(n)%AXIS_1d(:) )

             case('CX','CDX','CBFX')

                if ( ismaster ) then
                   GD1 = ( D1 - 2*xhalo ) * nprocs_x_out + 2*xhalo
                endif
                call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                call gather_x( ismaster,               &
                               bcast_,                 &
                               ainfo_all(n)%datatype,  &
                               nprocs_x_out,           &
                               nprocs_y_out,           &
                               D1, 0, xhalo, 0,        &
                               ainfo(n)%AXIS_1d(:),    &
                               ainfo_all(n)%AXIS_1d(:) )

             case('FX','FDX','FBFX')

                if ( ismaster ) then
                   GD1 = ( D1 - 2*xhalo - 1 ) * nprocs_x_out + 2*xhalo + 1
                endif
                call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                call gather_x( ismaster,               &
                               bcast_,                 &
                               ainfo_all(n)%datatype,  &
                               nprocs_x_out,           &
                               nprocs_y_out,           &
                               D1, 0, xhalo, 1,        &
                               ainfo(n)%AXIS_1d(:),    &
                               ainfo_all(n)%AXIS_1d(:) )

             case('CY','CDY','CBFY')

                if ( ismaster ) then
                   GD1 = ( D1 - 2*yhalo ) * nprocs_y_out + 2*yhalo
                endif
                call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                call gather_y( ismaster,               &
                               bcast_,                 &
                               ainfo_all(n)%datatype,  &
                               nprocs_x_out,           &
                               nprocs_y_out,           &
                               D1, 0, yhalo, 0,        &
                               ainfo(n)%AXIS_1d(:),    &
                               ainfo_all(n)%AXIS_1d(:) )

             case('FY','FDY','FBFY')

                if ( ismaster ) then
                   GD1 = ( D1 - 2*yhalo - 1 ) * nprocs_y_out + 2*yhalo + 1
                endif
                call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                call gather_y( ismaster,               &
                               bcast_,                 &
                               ainfo_all(n)%datatype,  &
                               nprocs_x_out,           &
                               nprocs_y_out,           &
                               D1, 0, yhalo, 1,        &
                               ainfo(n)%AXIS_1d(:),    &
                               ainfo_all(n)%AXIS_1d(:) )

             case default

                allocate( ainfo_all(n)%AXIS_1d( D1 ) )

                ainfo_all(n)%dim_size(1) = D1
                ainfo_all(n)%AXIS_1d(:) = ainfo(n)%AXIS_1d(:)

             endselect

          case( 2 )

             D1 = size( ainfo(n)%AXIS_2d(:,:), 1 )
             D2 = size( ainfo(n)%AXIS_2d(:,:), 2 )

             if ( index( ainfo_all(n)%varname, '_bnds' ) == 0 ) then

                select case( ainfo_all(n)%varname )
                case('lon_uy','lat_uy','cell_area_uy')

                   if ( ismaster ) then
                      GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                      GD2 = D2 * nprocs_y_out
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   call gather_xy( ismaster,                 &
                                   bcast_,                   &
                                   ainfo_all(n)%datatype,    &
                                   nprocs_x_out,             &
                                   nprocs_y_out,             &
                                   D1, 1, 0, 0,              &
                                   D2, 0, 0, 0,              &
                                   ainfo(n)%AXIS_2d(:,:),    &
                                   ainfo_all(n)%AXIS_2d(:,:) )

                case('lon_xv','lat_xv','cell_area_xv')

                   if ( ismaster ) then
                      GD1 = D1 * nprocs_x_out
                      GD2 = ( D2 - 1 ) * nprocs_y_out + 1
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   call gather_xy( ismaster,                 &
                                   bcast_,                   &
                                   ainfo_all(n)%datatype,    &
                                   nprocs_x_out,             &
                                   nprocs_y_out,             &
                                   D1, 0, 0, 0,              &
                                   D2, 1, 0, 0,              &
                                   ainfo(n)%AXIS_2d(:,:),    &
                                   ainfo_all(n)%AXIS_2d(:,:) )

                case('lon_uv','lat_uv')

                   if ( ismaster ) then
                      GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                      GD2 = ( D2 - 1 ) * nprocs_y_out + 1
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   call gather_xy( ismaster,                 &
                                   bcast_,                   &
                                   ainfo_all(n)%datatype,    &
                                   nprocs_x_out,             &
                                   nprocs_y_out,             &
                                   D1, 1, 0, 0,              &
                                   D2, 1, 0, 0,              &
                                   ainfo(n)%AXIS_2d(:,:),    &
                                   ainfo_all(n)%AXIS_2d(:,:) )

                case default

                   if ( ismaster ) then
                      GD1 = D1 * nprocs_x_out
                      GD2 = D2 * nprocs_y_out
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   call gather_xy( ismaster,                 &
                                   bcast_,                   &
                                   ainfo_all(n)%datatype,    &
                                   nprocs_x_out,             &
                                   nprocs_y_out,             &
                                   D1, 0, 0, 0,              &
                                   D2, 0, 0, 0,              &
                                   ainfo(n)%AXIS_2d(:,:),    &
                                   ainfo_all(n)%AXIS_2d(:,:) )

                endselect

             else

                select case( ainfo_all(n)%varname )
                case('x_bnds')

                   if ( ismaster ) then
                      GD1 = D1
                      GD2 = D2 * nprocs_x_out
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   call gather_nx( ismaster,                 &
                                   bcast_,                   &
                                   ainfo_all(n)%datatype,    &
                                   nprocs_x_out,             &
                                   nprocs_y_out,             &
                                   D1,                       &
                                   D2, 0, 0, 0,              &
                                   ainfo(n)%AXIS_2d(:,:),    &
                                   ainfo_all(n)%AXIS_2d(:,:) )

                case('xh_bnds')

                   if ( ismaster ) then
                      GD1 = D1
                      GD2 = ( D2 - 1 ) * nprocs_x_out + 1
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   call gather_nx( ismaster,                 &
                                   bcast_,                   &
                                   ainfo_all(n)%datatype,    &
                                   nprocs_x_out,             &
                                   nprocs_y_out,             &
                                   D1,                       &
                                   D2, 1, 0, 0,              &
                                   ainfo(n)%AXIS_2d(:,:),    &
                                   ainfo_all(n)%AXIS_2d(:,:) )

                case('y_bnds')

                   if ( ismaster ) then
                      GD1 = D1
                      GD2 = D2 * nprocs_y_out
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   call gather_ny( ismaster,                 &
                                   bcast_,                   &
                                   ainfo_all(n)%datatype,    &
                                   nprocs_x_out,             &
                                   nprocs_y_out,             &
                                   D1,                       &
                                   D2, 0, 0, 0,              &
                                   ainfo(n)%AXIS_2d(:,:),    &
                                   ainfo_all(n)%AXIS_2d(:,:) )

                case('yh_bnds')

                   if ( ismaster ) then
                      GD1 = D1
                      GD2 = ( D2 - 1 ) * nprocs_y_out + 1
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   call gather_ny( ismaster,                 &
                                   bcast_,                   &
                                   ainfo_all(n)%datatype,    &
                                   nprocs_x_out,             &
                                   nprocs_y_out,             &
                                   D1,                       &
                                   D2, 1, 0, 0,              &
                                   ainfo(n)%AXIS_2d(:,:),    &
                                   ainfo_all(n)%AXIS_2d(:,:) )

                case default

                   allocate( ainfo_all(n)%AXIS_2d( D1, D2 ) )

                   ainfo_all(n)%dim_size(1) = D1
                   ainfo_all(n)%dim_size(2) = D2
                   ainfo_all(n)%AXIS_2d(:,:) = ainfo(n)%AXIS_2d(:,:)

                endselect

             endif

          case( 3 )

             D1 = size( ainfo(n)%AXIS_3d(:,:,:), 1 )
             D2 = size( ainfo(n)%AXIS_3d(:,:,:), 2 )
             D3 = size( ainfo(n)%AXIS_3d(:,:,:), 3 )

             if ( ainfo_all(n)%transpose ) then

                select case( ainfo_all(n)%varname )
                case('height_uyz','height_uyw','cell_area_uyz_x','cell_area_uyw_x','cell_volume_uyz')

                   if ( ismaster ) then
                      GD1 = D1
                      GD2 = ( D2 - 1 ) * nprocs_x_out + 1
                      GD3 = D3 * nprocs_y_out
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD3, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, GD3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   call gather_nxy( ismaster,                   &
                                    bcast_,                     &
                                    ainfo_all(n)%datatype,      &
                                    nprocs_x_out,               &
                                    nprocs_y_out,               &
                                    D1,                         &
                                    D2, 1, 0, 0,                &
                                    D3, 0, 0, 0,                &
                                    ainfo(n)%AXIS_3d(:,:,:),    &
                                    ainfo_all(n)%AXIS_3d(:,:,:) )

                case('height_xvz','height_xvw','cell_area_xvz_y','cell_area_xvw_y','cell_volume_xvz')

                   if ( ismaster ) then
                      GD1 = D1
                      GD2 = D2 * nprocs_x_out
                      GD3 = ( D3 - 1 ) * nprocs_y_out + 1
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD3, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, GD3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   call gather_nxy( ismaster,                   &
                                    bcast_,                     &
                                    ainfo_all(n)%datatype,      &
                                    nprocs_x_out,               &
                                    nprocs_y_out,               &
                                    D1,                         &
                                    D2, 0, 0, 0,                &
                                    D3, 1, 0, 0,                &
                                    ainfo(n)%AXIS_3d(:,:,:),    &
                                    ainfo_all(n)%AXIS_3d(:,:,:) )

                case('height_uvz','height_uvw','cell_area_uvz_y','cell_area_uvz_x')

                   if ( ismaster ) then
                      GD1 = D1
                      GD2 = ( D2 - 1 ) * nprocs_x_out + 1
                      GD3 = ( D3 - 1 ) * nprocs_y_out + 1
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD3, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, GD3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   call gather_nxy( ismaster,                   &
                                    bcast_,                     &
                                    ainfo_all(n)%datatype,      &
                                    nprocs_x_out,               &
                                    nprocs_y_out,               &
                                    D1,                         &
                                    D2, 1, 0, 0,                &
                                    D3, 1, 0, 0,                &
                                    ainfo(n)%AXIS_3d(:,:,:),    &
                                    ainfo_all(n)%AXIS_3d(:,:,:) )

                case default

                   if ( ismaster ) then
                      GD1 = D1
                      GD2 = D2 * nprocs_x_out
                      GD3 = D3 * nprocs_y_out
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD3, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, GD3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   call gather_nxy( ismaster,                   &
                                    bcast_,                     &
                                    ainfo_all(n)%datatype,      &
                                    nprocs_x_out,               &
                                    nprocs_y_out,               &
                                    D1,                         &
                                    D2, 0, 0, 0,                &
                                    D3, 0, 0, 0,                &
                                    ainfo(n)%AXIS_3d(:,:,:),    &
                                    ainfo_all(n)%AXIS_3d(:,:,:) )

                endselect

             else

                select case( ainfo_all(n)%varname )
                case('height_uyz','height_uyw','cell_area_uyz_x','cell_area_uyw_x','cell_volume_uyz')

                   if ( ismaster ) then
                      GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                      GD2 = D2 * nprocs_y_out
                      GD3 = D3
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD3, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, D3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   call gather_xyn( ismaster,                   &
                                    bcast_,                     &
                                    ainfo_all(n)%datatype,      &
                                    nprocs_x_out,               &
                                    nprocs_y_out,               &
                                    D1, 1, 0, 0,                &
                                    D2, 0, 0, 0,                &
                                    D3,                         &
                                    ainfo(n)%AXIS_3d(:,:,:),    &
                                    ainfo_all(n)%AXIS_3d(:,:,:) )

                case('height_xvz','height_xvw','cell_area_xvz_y','cell_area_xvw_y','cell_volume_xvz')

                   if ( ismaster ) then
                      GD1 = D1 * nprocs_x_out
                      GD2 = ( D2 - 1 ) * nprocs_y_out + 1
                      GD3 = D3
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD3, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, D3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   call gather_xyn( ismaster,                   &
                                    bcast_,                     &
                                    ainfo_all(n)%datatype,      &
                                    nprocs_x_out,               &
                                    nprocs_y_out,               &
                                    D1, 0, 0, 0,                &
                                    D2, 1, 0, 0,                &
                                    D3,                         &
                                    ainfo(n)%AXIS_3d(:,:,:),    &
                                    ainfo_all(n)%AXIS_3d(:,:,:) )

                case('height_uvz','height_uvw','cell_area_uvz_y','cell_area_uvz_x')

                   if ( ismaster ) then
                      GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                      GD2 = ( D2 - 1 ) * nprocs_y_out + 1
                      GD3 = D3
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD3, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, D3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   call gather_xyn( ismaster,                   &
                                    bcast_,                     &
                                    ainfo_all(n)%datatype,      &
                                    nprocs_x_out,               &
                                    nprocs_y_out,               &
                                    D1, 1, 0, 0,                &
                                    D2, 1, 0, 0,                &
                                    D3,                         &
                                    ainfo(n)%AXIS_3d(:,:,:),    &
                                    ainfo_all(n)%AXIS_3d(:,:,:) )

                case default

                   if ( ismaster ) then
                      GD1 = D1 * nprocs_x_out
                      GD2 = D2 * nprocs_y_out
                      GD3 = D3
                   endif
                   call MPI_BCAST( GD1, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD2, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
                   call MPI_BCAST( GD3, 1, MPI_INTEGER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, D3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   call gather_xyn( ismaster,                   &
                                    bcast_,                     &
                                    ainfo_all(n)%datatype,      &
                                    nprocs_x_out,               &
                                    nprocs_y_out,               &
                                    D1, 0, 0, 0,                &
                                    D2, 0, 0, 0,                &
                                    D3,                         &
                                    ainfo(n)%AXIS_3d(:,:,:),    &
                                    ainfo_all(n)%AXIS_3d(:,:,:) )

                endselect

             endif

          endselect
       enddo

    else
       ! copy axis information
       ainfo_all(:) = ainfo(:)

    endif

    return
  end subroutine SNO_comm_globalaxis

  !-----------------------------------------------------------------------------
  subroutine SNO_comm_globalvars( &
       ismaster,      &
       output_single, &
       nprocs_x_out,  &
       nprocs_y_out,  &
       dinfo,         &
       dinfo_all,     &
       bcast          )
    use mod_sno_h, only: &
       iteminfo
    implicit none

    logical,          intent(in)    :: ismaster                              ! master rank process?
    logical,          intent(in)    :: output_single                         ! output single file when using MPI?
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (input)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (input)
    type(iteminfo),   intent(in)    :: dinfo                                 ! variable information               (input)
    type(iteminfo),   intent(out)   :: dinfo_all                             ! variable information               (output)

    logical,          intent(in), optional :: bcast                          ! use broadcast?

    integer  :: D1, D2, D3
    integer  :: GD1, GD2, GD3
    integer  :: t

    logical  :: bcast_
    !---------------------------------------------------------------------------

    if ( output_single ) then

       if ( present( bcast ) ) then
          bcast_ = .true.
       else
          bcast_ = .false.
       endif

       ! set variable information
       dinfo_all%varname       = dinfo%varname
       dinfo_all%description   = dinfo%description
       dinfo_all%units         = dinfo%units
       dinfo_all%standard_name = dinfo%standard_name
       dinfo_all%datatype      = dinfo%datatype
       dinfo_all%dim_rank      = dinfo%dim_rank
       dinfo_all%transpose     = dinfo%transpose
       dinfo_all%dim_name(:)   = dinfo%dim_name(:)
       dinfo_all%natts         = dinfo%natts
       dinfo_all%att_name(:)   = dinfo%att_name(:)
       dinfo_all%att_type(:)   = dinfo%att_type(:)
       dinfo_all%att_len (:)   = dinfo%att_len (:)
       dinfo_all%atts    (:)   = dinfo%atts    (:)
       dinfo_all%step_nmax     = dinfo%step_nmax
       do t = 1, dinfo%step_nmax
          dinfo_all%time_start(t) = dinfo%time_start(t)
          dinfo_all%time_end  (t) = dinfo%time_end  (t)
       enddo
       dinfo_all%dt            = dinfo%dt
       dinfo_all%time_units    = dinfo%time_units
       dinfo_all%calendar      = dinfo%calendar

       select case( dinfo_all%dim_rank )
       case( 1 )

          D1 = size( dinfo%VAR_1d(:), 1 )

          select case( dinfo_all%dim_name(1) )
          case('x')

             GD1 = D1 * nprocs_x_out

             allocate( dinfo_all%VAR_1d( GD1 ) )
             dinfo_all%dim_size(1) = GD1

             call gather_x( ismaster,           &
                            bcast_,             &
                            dinfo_all%datatype, &
                            nprocs_x_out,       &
                            nprocs_y_out,       &
                            D1, 0, 0, 0,        &
                            dinfo%VAR_1d(:),    &
                            dinfo_all%VAR_1d(:) )

          case('y')

             GD1 = D1 * nprocs_y_out

             allocate( dinfo_all%VAR_1d( GD1 ) )
             dinfo_all%dim_size(1) = GD1

             call gather_y( ismaster,           &
                            bcast_,             &
                            dinfo_all%datatype, &
                            nprocs_x_out,       &
                            nprocs_y_out,       &
                            D1, 0, 0, 0,        &
                            dinfo%VAR_1d(:),    &
                            dinfo_all%VAR_1d(:) )

          case default

             allocate( dinfo_all%VAR_1d( D1 ) )

             dinfo_all%dim_size(1) = D1
             dinfo_all%VAR_1d(:) = dinfo%VAR_1d(:)

          endselect

       case( 2 )

          D1 = size( dinfo%VAR_2d(:,:), 1 )
          D2 = size( dinfo%VAR_2d(:,:), 2 )

          if ( dinfo_all%dim_name(1) == 'lon' .AND. &
               dinfo_all%dim_name(2) == 'lat'       ) then

             GD1 = ( D1 - 1 ) * nprocs_x_out + 1
             GD2 = ( D2 - 1 ) * nprocs_y_out + 1

             allocate( dinfo_all%VAR_2d( GD1, GD2) )
             dinfo_all%dim_size(1) = GD1
             dinfo_all%dim_size(2) = GD2

             call gather_xy( ismaster,             &
                             bcast_,               &
                             dinfo_all%datatype,   &
                             nprocs_x_out,         &
                             nprocs_y_out,         &
                             D1, 1, 0, 0,          &
                             D2, 1, 0, 0,          &
                             dinfo%VAR_2d(:,:),    &
                             dinfo_all%VAR_2d(:,:) )

          else

             GD1 = D1 * nprocs_x_out
             GD2 = D2 * nprocs_y_out

             allocate( dinfo_all%VAR_2d( GD1, GD2) )
             dinfo_all%dim_size(1) = GD1
             dinfo_all%dim_size(2) = GD2

             call gather_xy( ismaster,             &
                             bcast_,               &
                             dinfo_all%datatype,   &
                             nprocs_x_out,         &
                             nprocs_y_out,         &
                             D1, 0, 0, 0,          &
                             D2, 0, 0, 0,          &
                             dinfo%VAR_2d(:,:),    &
                             dinfo_all%VAR_2d(:,:) )

          endif

       case( 3 )

          D1 = size( dinfo%VAR_3d(:,:,:), 1 )
          D2 = size( dinfo%VAR_3d(:,:,:), 2 )
          D3 = size( dinfo%VAR_3d(:,:,:), 3 )

          if ( dinfo_all%transpose ) then

             if ( dinfo_all%dim_name(1) == 'lon' .AND. &
                  dinfo_all%dim_name(2) == 'lat'       ) then

                GD1 = D1
                GD2 = ( D2 - 1 ) * nprocs_x_out + 1
                GD3 = ( D3 - 1 ) * nprocs_y_out + 1

                allocate( dinfo_all%VAR_3d( GD1, GD2, GD3 ) )
                dinfo_all%dim_size(1) = GD1
                dinfo_all%dim_size(2) = GD2
                dinfo_all%dim_size(3) = GD3

                call gather_nxy( ismaster,               &
                                 bcast_,                 &
                                 dinfo_all%datatype,     &
                                 nprocs_x_out,           &
                                 nprocs_y_out,           &
                                 D1,                     &
                                 D2, 1, 0, 0,            &
                                 D3, 1, 0, 0,            &
                                 dinfo%VAR_3d(:,:,:),    &
                                 dinfo_all%VAR_3d(:,:,:) )

             else

                GD1 = D1
                GD2 = D2 * nprocs_x_out
                GD3 = D3 * nprocs_y_out

                allocate( dinfo_all%VAR_3d( GD1, GD2, GD3 ) )
                dinfo_all%dim_size(1) = GD1
                dinfo_all%dim_size(2) = GD2
                dinfo_all%dim_size(3) = GD3

                call gather_nxy( ismaster,               &
                                 bcast_,                 &
                                 dinfo_all%datatype,     &
                                 nprocs_x_out,           &
                                 nprocs_y_out,           &
                                 D1,                     &
                                 D2, 0, 0, 0,            &
                                 D3, 0, 0, 0,            &
                                 dinfo%VAR_3d(:,:,:),    &
                                 dinfo_all%VAR_3d(:,:,:) )

             endif

          else

             if ( dinfo_all%dim_name(1) == 'lon' .AND. &
                  dinfo_all%dim_name(2) == 'lat'       ) then

                GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                GD2 = ( D2 - 1 ) * nprocs_y_out + 1
                GD3 = D3

                allocate( dinfo_all%VAR_3d( GD1, GD2, GD3 ) )
                dinfo_all%dim_size(1) = GD1
                dinfo_all%dim_size(2) = GD2
                dinfo_all%dim_size(3) = GD3

                call gather_xyn( ismaster,               &
                                 bcast_,                 &
                                 dinfo_all%datatype,     &
                                 nprocs_x_out,           &
                                 nprocs_y_out,           &
                                 D1, 1, 0, 0,            &
                                 D2, 1, 0, 0,            &
                                 D3,                     &
                                 dinfo%VAR_3d(:,:,:),    &
                                 dinfo_all%VAR_3d(:,:,:) )

             else

                GD1 = D1 * nprocs_x_out
                GD2 = D2 * nprocs_y_out
                GD3 = D3

                allocate( dinfo_all%VAR_3d( GD1, GD2, GD3 ) )
                dinfo_all%dim_size(1) = GD1
                dinfo_all%dim_size(2) = GD2
                dinfo_all%dim_size(3) = GD3

                call gather_xyn( ismaster,               &
                                 bcast_,                 &
                                 dinfo_all%datatype,     &
                                 nprocs_x_out,           &
                                 nprocs_y_out,           &
                                 D1, 0, 0, 0,            &
                                 D2, 0, 0, 0,            &
                                 D3,                     &
                                 dinfo%VAR_3d(:,:,:),    &
                                 dinfo_all%VAR_3d(:,:,:) )

             endif

          endif

       endselect

    else
       ! copy variable information
       dinfo_all = dinfo

    endif

    return
  end subroutine SNO_comm_globalvars

  !-----------------------------------------------------------------------------
  subroutine gather_x( &
       ismaster,       &
       bcast,          &
       FILE_RP,        &
       nprocs_x_out,   &
       nprocs_y_out,   &
       D1, H1, L1, M1, &
       din,            &
       dout            )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    logical,  intent(in)  :: bcast
    integer,  intent(in)  :: FILE_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1              ! local size of 1st dimension
    integer,  intent(in)  :: H1              ! subtraction of receive array count at D1 for half-level var. [0-1]
    integer,  intent(in)  :: L1              ! subtraction of output array count at D1 for var. with halo [halo size]
    integer,  intent(in)  :: M1              ! subtraction of outout array count at D1 for half-level var. with halo [0-1]

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(SP), allocatable :: recv_SP(:)
    real(DP), allocatable :: recv_DP(:)

    integer  :: sendpcnt
    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: iloc

    integer  :: i, m, p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D1
       else
          recvicnt(px) = D1 - H1
       endif

       recvpcnt(p) = recvicnt(px)

       if ( p == 1 ) then
          recvploc(p) = 0
       else
          recvploc(p) = sum( recvpcnt(1:p-1) )
       endif

       p = p + 1
    enddo
    enddo

    sendpcnt = D1

    select case( FILE_RP )
    case( FILE_REAL4 )
       allocate( recv_SP( sum( recvpcnt(:) ) ) )

       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHERV( send_SP(:),           &
                         sendpcnt,             &
                         MPI_REAL,             &
                         recv_SP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( FILE_REAL8 )
       allocate( recv_DP( sum( recvpcnt(:) ) ) )

       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHERV( send_DP(:),           &
                         sendpcnt,             &
                         MPI_DOUBLE_PRECISION, &
                         recv_DP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_DOUBLE_PRECISION, &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          if ( px == 1 ) then
             iloc = 0
          else
             iloc = sum( recvicnt(1:px-1) ) - ( 2 * L1 + M1 ) * ( px - 1 )
          endif

          select case( FILE_RP )
          case( FILE_REAL4 )
             ! rearrangement of data array
             m = 1
             do i = 1, recvicnt(px)
                dout( i + iloc ) = real( recv_SP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
          case( FILE_REAL8 )
             ! rearrangement of data array
             m = 1
             do i = 1, recvicnt(px)
                dout( i + iloc ) = real( recv_DP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
          endselect

          p = p + 1
       enddo
       enddo

    endif

    if ( bcast ) then

       sendpcnt = size( dout(:) )

       select case( RP )
       case( SP )
          call MPI_BCAST( dout(:),              &
                          sendpcnt,             &
                          MPI_REAL,             &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       case( DP )
          call MPI_BCAST( dout(:),              &
                          sendpcnt,             &
                          MPI_DOUBLE_PRECISION, &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       endselect

    endif

    return
  end subroutine gather_x

  !-----------------------------------------------------------------------------
  subroutine gather_y( &
       ismaster,       &
       bcast,          &
       FILE_RP,        &
       nprocs_x_out,   &
       nprocs_y_out,   &
       D1, H1, L1, M1, &
       din,            &
       dout            )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    logical,  intent(in)  :: bcast
    integer,  intent(in)  :: FILE_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1              ! local size of 1st dimension
    integer,  intent(in)  :: H1              ! subtraction of receive array count at D1 for half-level var. [0-1]
    integer,  intent(in)  :: L1              ! subtraction of output array count at D1 for var. with halo [halo size]
    integer,  intent(in)  :: M1              ! subtraction of outout array count at D1 for half-level var. with halo [0-1]

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(SP), allocatable :: recv_SP(:)
    real(DP), allocatable :: recv_DP(:)

    integer  :: sendpcnt
    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: jloc

    integer  :: j, m, p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( py == 1 ) then
          recvjcnt(py) = D1
       else
          recvjcnt(py) = D1 - H1
       endif

       recvpcnt(p) = recvjcnt(py)

       if ( p == 1 ) then
          recvploc(p) = 0
       else
          recvploc(p) = sum( recvpcnt(1:p-1) )
       endif

       p = p + 1
    enddo
    enddo

    sendpcnt = D1

    select case( FILE_RP )
    case( FILE_REAL4 )
       allocate( recv_SP( sum( recvpcnt(:) ) ) )

       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHERV( send_SP(:),           &
                         sendpcnt,             &
                         MPI_REAL,             &
                         recv_SP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( FILE_REAL8 )
       allocate( recv_DP( sum( recvpcnt(:) ) ) )

       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHERV( send_DP(:),           &
                         sendpcnt,             &
                         MPI_DOUBLE_PRECISION, &
                         recv_DP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_DOUBLE_PRECISION, &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          if ( py == 1 ) then
             jloc = 0
          else
             jloc = sum( recvjcnt(1:py-1) ) - ( 2 * L1 + M1 ) * ( py - 1 )
          endif

          select case( FILE_RP )
          case( FILE_REAL4 )
             ! rearrangement of data array
             m = 1
             do j = 1, recvjcnt(py)
                dout( j + jloc ) = real( recv_SP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
          case( FILE_REAL8 )
             ! rearrangement of data array
             m = 1
             do j = 1, recvjcnt(py)
                dout( j + jloc ) = real( recv_DP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
          endselect

          p = p + 1
       enddo
       enddo

    endif

    if ( bcast ) then

       sendpcnt = size( dout(:) )

       select case( RP )
       case( SP )
          call MPI_BCAST( dout(:),              &
                          sendpcnt,             &
                          MPI_REAL,             &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       case( DP )
          call MPI_BCAST( dout(:),              &
                          sendpcnt,             &
                          MPI_DOUBLE_PRECISION, &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       endselect

    endif

    return
  end subroutine gather_y

  !-----------------------------------------------------------------------------
  subroutine gather_xy( &
       ismaster,       &
       bcast,          &
       FILE_RP,        &
       nprocs_x_out,   &
       nprocs_y_out,   &
       D1, H1, L1, M1, &
       D2, H2, L2, M2, &
       din,            &
       dout            )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    logical,  intent(in)  :: bcast
    integer,  intent(in)  :: FILE_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1              ! local size of 1st dimension
    integer,  intent(in)  :: H1              ! subtraction of receive array count at D1 for half-level var. [0-1]
    integer,  intent(in)  :: L1              ! subtraction of output array count at D1 for var. with halo [halo size]
    integer,  intent(in)  :: M1              ! subtraction of outout array count at D1 for half-level var. with halo [0-1]
    integer,  intent(in)  :: D2              ! local size of 2nd dimension
    integer,  intent(in)  :: H2              ! subtraction of receive array count at D2 for half-level var. [0-1]
    integer,  intent(in)  :: L2              ! subtraction of output array count at D2 for var. with halo [halo size]
    integer,  intent(in)  :: M2              ! subtraction of outout array count at D2 for half-level var. with halo [0-1]

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(SP), allocatable :: recv_SP(:)
    real(DP), allocatable :: recv_DP(:)

    integer  :: sendpcnt
    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: iloc, jloc

    integer  :: i, j, m, p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D1
       else
          recvicnt(px) = D1 - H1
       endif
       if ( py == 1 ) then
          recvjcnt(py) = D2
       else
          recvjcnt(py) = D2 - H2
       endif

       recvpcnt(p) = recvicnt(px) * recvjcnt(py)

       if ( p == 1 ) then
          recvploc(p) = 0
       else
          recvploc(p) = sum( recvpcnt(1:p-1) )
       endif

       p = p + 1
    enddo
    enddo

    sendpcnt = D1 * D2

    select case( FILE_RP )
    case( FILE_REAL4 )
       allocate( recv_SP( sum( recvpcnt(:) ) ) )

       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:),         &
                         sendpcnt,             &
                         MPI_REAL,             &
                         recv_SP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( FILE_REAL8 )
       allocate( recv_DP( sum( recvpcnt(:) ) ) )

       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:),         &
                         sendpcnt,             &
                         MPI_DOUBLE_PRECISION, &
                         recv_DP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_DOUBLE_PRECISION, &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          if ( px == 1 ) then
             iloc = 0
          else
             iloc = sum( recvicnt(1:px-1) ) - ( 2 * L1 + M1 ) * ( px - 1 )
          endif
          if ( py == 1 ) then
             jloc = 0
          else
             jloc = sum( recvjcnt(1:py-1) ) - ( 2 * L2 + M2 ) * ( py - 1 )
          endif

          select case( FILE_RP )
          case( FILE_REAL4 )
             ! rearrangement of data array
             m = 1
             do j = 1, recvjcnt(py)
             do i = 1, recvicnt(px)
                dout( i + iloc, &
                      j + jloc  ) = real( recv_SP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
          case( FILE_REAL8 )
             ! rearrangement of data array
             m = 1
             do j = 1, recvjcnt(py)
             do i = 1, recvicnt(px)
                dout( i + iloc, &
                      j + jloc  ) = real( recv_DP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
          endselect

          p = p + 1
       enddo
       enddo

    endif

    if ( bcast ) then

       sendpcnt = size( dout(:,:) )

       select case( RP )
       case( SP )
          call MPI_BCAST( dout(:,:),            &
                          sendpcnt,             &
                          MPI_REAL,             &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       case( DP )
          call MPI_BCAST( dout(:,:),            &
                          sendpcnt,             &
                          MPI_DOUBLE_PRECISION, &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       endselect

    endif

    return
  end subroutine gather_xy

  !-----------------------------------------------------------------------------
  subroutine gather_nx( &
       ismaster,       &
       bcast,          &
       FILE_RP,        &
       nprocs_x_out,   &
       nprocs_y_out,   &
       D1,             &
       D2, H2, L2, M2, &
       din,            &
       dout            )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    logical,  intent(in)  :: bcast
    integer,  intent(in)  :: FILE_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1              ! local size of 1st dimension
    integer,  intent(in)  :: D2              ! local size of 2nd dimension
    integer,  intent(in)  :: H2              ! subtraction of receive array count at D2 for half-level var. [0-1]
    integer,  intent(in)  :: L2              ! subtraction of output array count at D2 for var. with halo [halo size]
    integer,  intent(in)  :: M2              ! subtraction of outout array count at D2 for half-level var. with halo [0-1]

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(SP), allocatable :: recv_SP(:)
    real(DP), allocatable :: recv_DP(:)

    integer  :: sendpcnt
    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: iloc

    integer  :: v, i, m, p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D2
       else
          recvicnt(px) = D2 - H2
       endif

       recvpcnt(p) = D1 * recvicnt(px)

       if ( p == 1 ) then
          recvploc(p) = 0
       else
          recvploc(p) = sum( recvpcnt(1:p-1) )
       endif

       p = p + 1
    enddo
    enddo

    sendpcnt = D1 * D2

    select case( FILE_RP )
    case( FILE_REAL4 )
       allocate( recv_SP( sum( recvpcnt(:) ) ) )

       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:),         &
                         sendpcnt,             &
                         MPI_REAL,             &
                         recv_SP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( FILE_REAL8 )
       allocate( recv_DP( sum( recvpcnt(:) ) ) )

       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:),         &
                         sendpcnt,             &
                         MPI_DOUBLE_PRECISION, &
                         recv_DP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_DOUBLE_PRECISION, &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          if ( px == 1 ) then
             iloc = 0
          else
             iloc = sum( recvicnt(1:px-1) ) - ( 2 * L2 + M2 ) * ( px - 1 )
          endif

          select case( FILE_RP )
          case( FILE_REAL4 )
             ! rearrangement of data array
             m = 1
             do i = 1, recvicnt(px)
             do v = 1, D1
                dout( v,       &
                      i + iloc ) = real( recv_SP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
          case( FILE_REAL8 )
             ! rearrangement of data array
             m = 1
             do i = 1, recvicnt(px)
             do v = 1, D1
                dout( v,       &
                      i + iloc ) = real( recv_DP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
          endselect

          p = p + 1
       enddo
       enddo

    endif

    if ( bcast ) then

       sendpcnt = size( dout(:,:) )

       select case( RP )
       case( SP )
          call MPI_BCAST( dout(:,:),            &
                          sendpcnt,             &
                          MPI_REAL,             &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       case( DP )
          call MPI_BCAST( dout(:,:),            &
                          sendpcnt,             &
                          MPI_DOUBLE_PRECISION, &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       endselect

    endif

    return
  end subroutine gather_nx

  !-----------------------------------------------------------------------------
  subroutine gather_ny( &
       ismaster,       &
       bcast,          &
       FILE_RP,        &
       nprocs_x_out,   &
       nprocs_y_out,   &
       D1,             &
       D2, H2, L2, M2, &
       din,            &
       dout            )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    logical,  intent(in)  :: bcast
    integer,  intent(in)  :: FILE_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1              ! local size of 1st dimension
    integer,  intent(in)  :: D2              ! local size of 2nd dimension
    integer,  intent(in)  :: H2              ! subtraction of receive array count at D2 for half-level var. [0-1]
    integer,  intent(in)  :: L2              ! subtraction of output array count at D2 for var. with halo [halo size]
    integer,  intent(in)  :: M2              ! subtraction of outout array count at D2 for half-level var. with halo [0-1]

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(SP), allocatable :: recv_SP(:)
    real(DP), allocatable :: recv_DP(:)

    integer  :: sendpcnt
    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: jloc

    integer  :: v, j, m, p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( py == 1 ) then
          recvjcnt(py) = D2
       else
          recvjcnt(py) = D2 - H2
       endif

       recvpcnt(p) = D1 * recvjcnt(py)

       if ( p == 1 ) then
          recvploc(p) = 0
       else
          recvploc(p) = sum( recvpcnt(1:p-1) )
       endif

       p = p + 1
    enddo
    enddo

    sendpcnt = D1 * D2

    select case( FILE_RP )
    case( FILE_REAL4 )
       allocate( recv_SP( sum( recvpcnt(:) ) ) )

       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:),         &
                         sendpcnt,             &
                         MPI_REAL,             &
                         recv_SP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( FILE_REAL8 )
       allocate( recv_DP( sum( recvpcnt(:) ) ) )

       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:),         &
                         sendpcnt,             &
                         MPI_DOUBLE_PRECISION, &
                         recv_DP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_DOUBLE_PRECISION, &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          if ( py == 1 ) then
             jloc = 0
          else
             jloc = sum( recvjcnt(1:py-1) ) - ( 2 * L2 + M2 ) * ( py - 1 )
          endif

          select case( FILE_RP )
          case( FILE_REAL4 )
             ! rearrangement of data array
             m = 1
             do j = 1, recvjcnt(py)
             do v = 1, D1
                dout( v,       &
                      j + jloc ) = real( recv_SP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
          case( FILE_REAL8 )
             ! rearrangement of data array
             m = 1
             do j = 1, recvjcnt(py)
             do v = 1, D1
                dout( v,       &
                      j + jloc ) = real( recv_DP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
          endselect

          p = p + 1
       enddo
       enddo

    endif

    if ( bcast ) then

       sendpcnt = size( dout(:,:) )

       select case( RP )
       case( SP )
          call MPI_BCAST( dout(:,:),            &
                          sendpcnt,             &
                          MPI_REAL,             &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       case( DP )
          call MPI_BCAST( dout(:,:),            &
                          sendpcnt,             &
                          MPI_DOUBLE_PRECISION, &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       endselect

    endif

    return
  end subroutine gather_ny

  !-----------------------------------------------------------------------------
  subroutine gather_nxy( &
       ismaster,       &
       bcast,          &
       FILE_RP,        &
       nprocs_x_out,   &
       nprocs_y_out,   &
       D1,             &
       D2, H2, L2, M2, &
       D3, H3, L3, M3, &
       din,            &
       dout            )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    logical,  intent(in)  :: bcast
    integer,  intent(in)  :: FILE_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1              ! local size of 1st dimension
    integer,  intent(in)  :: D2              ! local size of 2nd dimension
    integer,  intent(in)  :: H2              ! subtraction of receive array count at D2 for half-level var. [0-1]
    integer,  intent(in)  :: L2              ! subtraction of output array count at D2 for var. with halo [halo size]
    integer,  intent(in)  :: M2              ! subtraction of outout array count at D2 for half-level var. with halo [0-1]
    integer,  intent(in)  :: D3              ! local size of 3rd dimension
    integer,  intent(in)  :: H3              ! subtraction of receive array count at D3 for half-level var. [0-1]
    integer,  intent(in)  :: L3              ! subtraction of output array count at D3 for var. with halo [halo size]
    integer,  intent(in)  :: M3              ! subtraction of outout array count at D3 for half-level var. with halo [0-1]

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(SP), allocatable :: recv_SP(:)
    real(DP), allocatable :: recv_DP(:)

    integer  :: sendpcnt
    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: iloc, jloc

    integer  :: i, j, k, m, p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D2
       else
          recvicnt(px) = D2 - H2
       endif
       if ( py == 1 ) then
          recvjcnt(py) = D3
       else
          recvjcnt(py) = D3 - H3
       endif

       recvpcnt(p) = D1 * recvicnt(px) * recvjcnt(py)

       if ( p == 1 ) then
          recvploc(p) = 0
       else
          recvploc(p) = sum( recvpcnt(1:p-1) )
       endif

       p = p + 1
    enddo
    enddo

    sendpcnt = D1 * D2 * D3

    select case( FILE_RP )
    case( FILE_REAL4 )
       allocate( recv_SP( sum( recvpcnt(:) ) ) )

       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:,:),       &
                         sendpcnt,             &
                         MPI_REAL,             &
                         recv_SP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( FILE_REAL8 )
       allocate( recv_DP( sum( recvpcnt(:) ) ) )

       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:,:),       &
                         sendpcnt,             &
                         MPI_DOUBLE_PRECISION, &
                         recv_DP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_DOUBLE_PRECISION, &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          if ( px == 1 ) then
             iloc = 0
          else
             iloc = sum( recvicnt(1:px-1) ) - ( 2 * L2 + M2 ) * ( px - 1 )
          endif
          if ( py == 1 ) then
             jloc = 0
          else
             jloc = sum( recvjcnt(1:py-1) ) - ( 2 * L3 + M3 ) * ( py - 1 )
          endif

          select case( FILE_RP )
          case( FILE_REAL4 )
             ! rearrangement of data array
             m = 1
             do j = 1, recvjcnt(py)
             do i = 1, recvicnt(px)
             do k = 1, D1
                dout( k,        &
                      i + iloc, &
                      j + jloc  ) = real( recv_SP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
             enddo
          case( FILE_REAL8 )
             ! rearrangement of data array
             m = 1
             do j = 1, recvjcnt(py)
             do i = 1, recvicnt(px)
             do k = 1, D1
                dout( k,        &
                      i + iloc, &
                      j + jloc  ) = real( recv_DP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
             enddo
          endselect

          p = p + 1
       enddo
       enddo

    endif

    if ( bcast ) then

       sendpcnt = size( dout(:,:,:) )

       select case( RP )
       case( SP )
          call MPI_BCAST( dout(:,:,:),          &
                          sendpcnt,             &
                          MPI_REAL,             &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       case( DP )
          call MPI_BCAST( dout(:,:,:),          &
                          sendpcnt,             &
                          MPI_DOUBLE_PRECISION, &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       endselect

    endif

    return
  end subroutine gather_nxy

  !-----------------------------------------------------------------------------
  subroutine gather_xyn( &
       ismaster,       &
       bcast,          &
       FILE_RP,        &
       nprocs_x_out,   &
       nprocs_y_out,   &
       D1, H1, L1, M1, &
       D2, H2, L2, M2, &
       D3,             &
       din,            &
       dout            )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    logical,  intent(in)  :: bcast
    integer,  intent(in)  :: FILE_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1              ! local size of 1st dimension
    integer,  intent(in)  :: H1              ! subtraction of receive array count at D1 for half-level var. [0-1]
    integer,  intent(in)  :: L1              ! subtraction of output array count at D1 for var. with halo [halo size]
    integer,  intent(in)  :: M1              ! subtraction of outout array count at D1 for half-level var. with halo [0-1]
    integer,  intent(in)  :: D2              ! local size of 2nd dimension
    integer,  intent(in)  :: H2              ! subtraction of receive array count at D2 for half-level var. [0-1]
    integer,  intent(in)  :: L2              ! subtraction of output array count at D2 for var. with halo [halo size]
    integer,  intent(in)  :: M2              ! subtraction of outout array count at D2 for half-level var. with halo [0-1]
    integer,  intent(in)  :: D3              ! local size of 3rd dimension

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(SP), allocatable :: recv_SP(:)
    real(DP), allocatable :: recv_DP(:)

    integer  :: sendpcnt
    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: iloc, jloc

    integer  :: i, j, k, m, p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D1
       else
          recvicnt(px) = D1 - H1
       endif
       if ( py == 1 ) then
          recvjcnt(py) = D2
       else
          recvjcnt(py) = D2 - H2
       endif

       recvpcnt(p) = recvicnt(px) * recvjcnt(py) * D3

       if ( p == 1 ) then
          recvploc(p) = 0
       else
          recvploc(p) = sum( recvpcnt(1:p-1) )
       endif

       p = p + 1
    enddo
    enddo

    sendpcnt = D1 * D2 * D3

    select case( FILE_RP )
    case( FILE_REAL4 )
       allocate( recv_SP( sum( recvpcnt(:) ) ) )

       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:,:),       &
                         sendpcnt,             &
                         MPI_REAL,             &
                         recv_SP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( FILE_REAL8 )
       allocate( recv_DP( sum( recvpcnt(:) ) ) )

       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:,:),       &
                         sendpcnt,             &
                         MPI_DOUBLE_PRECISION, &
                         recv_DP (:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_DOUBLE_PRECISION, &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          if ( px == 1 ) then
             iloc = 0
          else
             iloc = sum( recvicnt(1:px-1) ) - ( 2 * L1 + M1 ) * ( px - 1 )
          endif
          if ( py == 1 ) then
             jloc = 0
          else
             jloc = sum( recvjcnt(1:py-1) ) - ( 2 * L2 + M2 ) * ( py - 1 )
          endif

          select case( FILE_RP )
          case( FILE_REAL4 )
             ! rearrangement of data array
             m = 1
             do k = 1, D3
             do j = 1, recvjcnt(py)
             do i = 1, recvicnt(px)
                dout( i + iloc, &
                      j + jloc, &
                      k         ) = real( recv_SP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
             enddo
          case( FILE_REAL8 )
             ! rearrangement of data array
             m = 1
             do k = 1, D3
             do j = 1, recvjcnt(py)
             do i = 1, recvicnt(px)
                dout( i + iloc, &
                      j + jloc, &
                      k         ) = real( recv_DP( m + recvploc(p) ), kind=RP )
                m = m + 1
             enddo
             enddo
             enddo
          endselect

          p = p + 1
       enddo
       enddo

    endif

    if ( bcast ) then

       sendpcnt = size( dout(:,:,:) )

       select case( RP )
       case( SP )
          call MPI_BCAST( dout(:,:,:),          &
                          sendpcnt,             &
                          MPI_REAL,             &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       case( DP )
          call MPI_BCAST( dout(:,:,:),          &
                          sendpcnt,             &
                          MPI_DOUBLE_PRECISION, &
                          PRC_masterrank,       &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )
       endselect

    endif

    return
  end subroutine gather_xyn

end module mod_sno_comm
