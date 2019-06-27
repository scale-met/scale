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
  private :: gather_u
  private :: gather_y
  private :: gather_v

  private :: gather_x_halo
  private :: gather_u_halo
  private :: gather_y_halo
  private :: gather_v_halo

  ! 2D
  private :: gather_xy
  private :: gather_uy
  private :: gather_xv
  private :: gather_uv

  private :: gather_nx
  private :: gather_nu
  private :: gather_ny
  private :: gather_nv

  ! 3D
  private :: gather_zxy
  private :: gather_zuy
  private :: gather_zxv
  private :: gather_zuv

  private :: gather_xyz
  private :: gather_uyz
  private :: gather_xvz
  private :: gather_uvz

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
       ainfo_all      )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
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

    integer  :: D1, D2, D3
    integer  :: GD1, GD2, GD3
    integer  :: xhalo, yhalo
    integer  :: n
    !---------------------------------------------------------------------------

    if ( output_single ) then

       xhalo = hinfo%halosize(2)
       yhalo = hinfo%halosize(3)

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
             case('x','lon')

                GD1 = D1 * nprocs_x_out

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )
                   call gather_x( ismaster,               &
                                  SP,                     &
                                  nprocs_x_out,           &
                                  nprocs_y_out,           &
                                  D1,                     &
                                  ainfo(n)%AXIS_1d(:),    &
                                  ainfo_all(n)%AXIS_1d(:) )
                case( FILE_REAL8 )
                   call gather_x( ismaster,               &
                                  DP,                     &
                                  nprocs_x_out,           &
                                  nprocs_y_out,           &
                                  D1,                     &
                                  ainfo(n)%AXIS_1d(:),    &
                                  ainfo_all(n)%AXIS_1d(:) )
                endselect

             case('xh')

                GD1 = ( D1 - 1 ) * nprocs_x_out + 1

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )
                   call gather_u( ismaster,               &
                                  SP,                     &
                                  nprocs_x_out,           &
                                  nprocs_y_out,           &
                                  D1,                     &
                                  ainfo(n)%AXIS_1d(:),    &
                                  ainfo_all(n)%AXIS_1d(:) )
                case( FILE_REAL8 )
                   call gather_u( ismaster,               &
                                  DP,                     &
                                  nprocs_x_out,           &
                                  nprocs_y_out,           &
                                  D1,                     &
                                  ainfo(n)%AXIS_1d(:),    &
                                  ainfo_all(n)%AXIS_1d(:) )
                endselect

             case('y','lat')

                GD1 = D1 * nprocs_y_out

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )
                   call gather_y( ismaster,               &
                                  SP,                     &
                                  nprocs_x_out,           &
                                  nprocs_y_out,           &
                                  D1,                     &
                                  ainfo(n)%AXIS_1d(:),    &
                                  ainfo_all(n)%AXIS_1d(:) )
                case( FILE_REAL8 )
                   call gather_y( ismaster,               &
                                  DP,                     &
                                  nprocs_x_out,           &
                                  nprocs_y_out,           &
                                  D1,                     &
                                  ainfo(n)%AXIS_1d(:),    &
                                  ainfo_all(n)%AXIS_1d(:) )
                endselect

             case('yh')

                GD1 = ( D1 - 1 ) * nprocs_y_out + 1

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )
                   call gather_v( ismaster,               &
                                  SP,                     &
                                  nprocs_x_out,           &
                                  nprocs_y_out,           &
                                  D1,                     &
                                  ainfo(n)%AXIS_1d(:),    &
                                  ainfo_all(n)%AXIS_1d(:) )
                case( FILE_REAL8 )
                   call gather_v( ismaster,               &
                                  DP,                     &
                                  nprocs_x_out,           &
                                  nprocs_y_out,           &
                                  D1,                     &
                                  ainfo(n)%AXIS_1d(:),    &
                                  ainfo_all(n)%AXIS_1d(:) )
                endselect

             case('CX','CDX','CBFX')

                GD1 = ( D1 - 2*xhalo ) * nprocs_x_out + 2*xhalo

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )
                   call gather_x_halo( ismaster,               &
                                       SP,                     &
                                       nprocs_x_out,           &
                                       nprocs_y_out,           &
                                       xhalo,                  &
                                       D1,                     &
                                       ainfo(n)%AXIS_1d(:),    &
                                       ainfo_all(n)%AXIS_1d(:) )
                case( FILE_REAL8 )
                   call gather_x_halo( ismaster,               &
                                       DP,                     &
                                       nprocs_x_out,           &
                                       nprocs_y_out,           &
                                       xhalo,                  &
                                       D1,                     &
                                       ainfo(n)%AXIS_1d(:),    &
                                       ainfo_all(n)%AXIS_1d(:) )
                endselect

             case('FX','FDX','FBFX')

                GD1 = ( D1 - 2*xhalo - 1 ) * nprocs_x_out + 2*xhalo + 1

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )
                   call gather_u_halo( ismaster,               &
                                       SP,                     &
                                       nprocs_x_out,           &
                                       nprocs_y_out,           &
                                       xhalo,                  &
                                       D1,                     &
                                       ainfo(n)%AXIS_1d(:),    &
                                       ainfo_all(n)%AXIS_1d(:) )
                case( FILE_REAL8 )
                   call gather_u_halo( ismaster,               &
                                       DP,                     &
                                       nprocs_x_out,           &
                                       nprocs_y_out,           &
                                       xhalo,                  &
                                       D1,                     &
                                       ainfo(n)%AXIS_1d(:),    &
                                       ainfo_all(n)%AXIS_1d(:) )
                endselect

             case('CY','CDY','CBFY')

                GD1 = ( D1 - 2*yhalo ) * nprocs_y_out + 2*yhalo

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )
                   call gather_y_halo( ismaster,               &
                                       SP,                     &
                                       nprocs_x_out,           &
                                       nprocs_y_out,           &
                                       yhalo,                  &
                                       D1,                     &
                                       ainfo(n)%AXIS_1d(:),    &
                                       ainfo_all(n)%AXIS_1d(:) )
                case( FILE_REAL8 )
                   call gather_y_halo( ismaster,               &
                                       DP,                     &
                                       nprocs_x_out,           &
                                       nprocs_y_out,           &
                                       yhalo,                  &
                                       D1,                     &
                                       ainfo(n)%AXIS_1d(:),    &
                                       ainfo_all(n)%AXIS_1d(:) )
                endselect

             case('FY','FDY','FBFY')

                GD1 = ( D1 - 2*yhalo - 1 ) * nprocs_y_out + 2*yhalo + 1

                allocate( ainfo_all(n)%AXIS_1d( GD1 ) )
                ainfo_all(n)%dim_size(1) = GD1

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )
                   call gather_v_halo( ismaster,               &
                                       SP,                     &
                                       nprocs_x_out,           &
                                       nprocs_y_out,           &
                                       yhalo,                  &
                                       D1,                     &
                                       ainfo(n)%AXIS_1d(:),    &
                                       ainfo_all(n)%AXIS_1d(:) )
                case( FILE_REAL8 )
                   call gather_v_halo( ismaster,               &
                                       DP,                     &
                                       nprocs_x_out,           &
                                       nprocs_y_out,           &
                                       yhalo,                  &
                                       D1,                     &
                                       ainfo(n)%AXIS_1d(:),    &
                                       ainfo_all(n)%AXIS_1d(:) )
                endselect

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

                   GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                   GD2 = D2 * nprocs_y_out

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_uy( ismaster,                 &
                                      SP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   case( FILE_REAL8 )
                      call gather_uy( ismaster,                 &
                                      DP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   endselect

                case('lon_xv','lat_xv','cell_area_xv')

                   GD1 = D1 * nprocs_x_out
                   GD2 = ( D2 - 1 ) * nprocs_y_out + 1

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_xv( ismaster,                 &
                                      SP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   case( FILE_REAL8 )
                      call gather_xv( ismaster,                 &
                                      DP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   endselect

                case('lon_uv','lat_uv')

                   GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                   GD2 = ( D2 - 1 ) * nprocs_y_out + 1

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_uv( ismaster,                 &
                                      SP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   case( FILE_REAL8 )
                      call gather_uv( ismaster,                 &
                                      DP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   endselect

                case default

                   GD1 = D1 * nprocs_x_out
                   GD2 = D2 * nprocs_y_out

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_xy( ismaster,                 &
                                      SP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   case( FILE_REAL8 )
                      call gather_xy( ismaster,                 &
                                      DP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   endselect

                endselect

             else

                select case( ainfo_all(n)%varname )
                case('x_bnds')

                   GD1 = D1
                   GD2 = D2 * nprocs_x_out

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_nx( ismaster,                 &
                                      SP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   case( FILE_REAL8 )
                      call gather_nx( ismaster,                 &
                                      DP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   endselect

                case('xh_bnds')

                   GD1 = D1
                   GD2 = ( D2 - 1 ) * nprocs_x_out + 1

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_nu( ismaster,                 &
                                      SP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   case( FILE_REAL8 )
                      call gather_nu( ismaster,                 &
                                      DP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   endselect

                case('y_bnds')

                   GD1 = D1
                   GD2 = D2 * nprocs_y_out

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_ny( ismaster,                 &
                                      SP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   case( FILE_REAL8 )
                      call gather_ny( ismaster,                 &
                                      DP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   endselect

                case('yh_bnds')

                   GD1 = D1
                   GD2 = ( D2 - 1 ) * nprocs_y_out + 1

                   allocate( ainfo_all(n)%AXIS_2d( GD1, GD2 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_nv( ismaster,                 &
                                      SP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   case( FILE_REAL8 )
                      call gather_nv( ismaster,                 &
                                      DP,                       &
                                      nprocs_x_out,             &
                                      nprocs_y_out,             &
                                      D1, D2,                   &
                                      ainfo(n)%AXIS_2d(:,:),    &
                                      ainfo_all(n)%AXIS_2d(:,:) )
                   endselect

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

                   GD1 = D1
                   GD2 = ( D2 - 1 ) * nprocs_x_out + 1
                   GD3 = D3 * nprocs_y_out

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, GD3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_zuy( ismaster,                   &
                                       SP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   case( FILE_REAL8 )
                      call gather_zuy( ismaster,                   &
                                       DP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   endselect

                case('height_xvz','height_xvw','cell_area_xvz_y','cell_area_xvw_y','cell_volume_xvz')

                   GD1 = D1
                   GD2 = D2 * nprocs_x_out
                   GD3 = ( D3 - 1 ) * nprocs_y_out + 1

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, GD3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_zxv( ismaster,                   &
                                       SP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   case( FILE_REAL8 )
                      call gather_zxv( ismaster,                   &
                                       DP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   endselect

                case('height_uvz','height_uvw','cell_area_uvz_y','cell_area_uvz_x')

                   GD1 = D1
                   GD2 = ( D2 - 1 ) * nprocs_x_out + 1
                   GD3 = ( D3 - 1 ) * nprocs_y_out + 1

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, GD3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_zuv( ismaster,                   &
                                       SP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   case( FILE_REAL8 )
                      call gather_zuv( ismaster,                   &
                                       DP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   endselect

                case default

                   GD1 = D1
                   GD2 = D2 * nprocs_x_out
                   GD3 = D3 * nprocs_y_out

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, GD3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_zxy( ismaster,                   &
                                       SP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   case( FILE_REAL8 )
                      call gather_zxy( ismaster,                   &
                                       DP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   endselect

                endselect

             else

                select case( ainfo_all(n)%varname )
                case('height_uyz','height_uyw','cell_area_uyz_x','cell_area_uyw_x','cell_volume_uyz')

                   GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                   GD2 = D2 * nprocs_y_out
                   GD3 = D3

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, D3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_uyz( ismaster,                   &
                                       SP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   case( FILE_REAL8 )
                      call gather_uyz( ismaster,                   &
                                       DP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   endselect

                case('height_xvz','height_xvw','cell_area_xvz_y','cell_area_xvw_y','cell_volume_xvz')

                   GD1 = D1 * nprocs_x_out
                   GD2 = ( D2 - 1 ) * nprocs_y_out + 1
                   GD3 = D3

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, D3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_xvz( ismaster,                   &
                                       SP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   case( FILE_REAL8 )
                      call gather_xvz( ismaster,                   &
                                       DP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   endselect

                case('height_uvz','height_uvw','cell_area_uvz_y','cell_area_uvz_x')

                   GD1 = ( D1 - 1 ) * nprocs_x_out + 1
                   GD2 = ( D2 - 1 ) * nprocs_y_out + 1
                   GD3 = D3

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, D3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_uvz( ismaster,                   &
                                       SP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   case( FILE_REAL8 )
                      call gather_uvz( ismaster,                   &
                                       DP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   endselect

                case default

                   GD1 = D1 * nprocs_x_out
                   GD2 = D2 * nprocs_y_out
                   GD3 = D3

                   allocate( ainfo_all(n)%AXIS_3d( GD1, GD2, D3 ) )
                   ainfo_all(n)%dim_size(1) = GD1
                   ainfo_all(n)%dim_size(2) = GD2
                   ainfo_all(n)%dim_size(3) = GD3

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )
                      call gather_xyz( ismaster,                   &
                                       SP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   case( FILE_REAL8 )
                      call gather_xyz( ismaster,                   &
                                       DP,                         &
                                       nprocs_x_out,               &
                                       nprocs_y_out,               &
                                       D1, D2, D3,                 &
                                       ainfo(n)%AXIS_3d(:,:,:),    &
                                       ainfo_all(n)%AXIS_3d(:,:,:) )
                   endselect

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
       dinfo_all      )
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use mod_sno_h, only: &
       iteminfo
    implicit none

    logical,          intent(in)    :: ismaster                              ! master rank process?
    logical,          intent(in)    :: output_single                         ! output single file when using MPI?
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (input)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (input)
    type(iteminfo),   intent(in)    :: dinfo                                 ! variable information               (input)
    type(iteminfo),   intent(out)   :: dinfo_all                             ! variable information               (output)

    integer  :: D1, D2, D3
    integer  :: GD1, GD2, GD3
    integer  :: t
    !---------------------------------------------------------------------------

    if ( output_single ) then

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
          case('x','lon')

             GD1 = D1 * nprocs_x_out

             allocate( dinfo_all%VAR_1d( GD1 ) )
             dinfo_all%dim_size(1) = GD1

             select case( dinfo_all%datatype )
             case( FILE_REAL4 )
                call gather_x( ismaster,           &
                               SP,                 &
                               nprocs_x_out,       &
                               nprocs_y_out,       &
                               D1,                 &
                               dinfo%VAR_1d(:),    &
                               dinfo_all%VAR_1d(:) )
             case( FILE_REAL8 )
                call gather_x( ismaster,           &
                               DP,                 &
                               nprocs_x_out,       &
                               nprocs_y_out,       &
                               D1,                 &
                               dinfo%VAR_1d(:),    &
                               dinfo_all%VAR_1d(:) )
             endselect

          case('y','lat')

             GD1 = D1 * nprocs_y_out

             allocate( dinfo_all%VAR_1d( GD1 ) )
             dinfo_all%dim_size(1) = GD1

             select case( dinfo_all%datatype )
             case( FILE_REAL4 )
                call gather_y( ismaster,           &
                               SP,                 &
                               nprocs_x_out,       &
                               nprocs_y_out,       &
                               D1,                 &
                               dinfo%VAR_1d(:),    &
                               dinfo_all%VAR_1d(:) )
             case( FILE_REAL8 )
                call gather_y( ismaster,           &
                               DP,                 &
                               nprocs_x_out,       &
                               nprocs_y_out,       &
                               D1,                 &
                               dinfo%VAR_1d(:),    &
                               dinfo_all%VAR_1d(:) )
             endselect

          case default

             allocate( dinfo_all%VAR_1d( D1 ) )

             dinfo_all%dim_size(1) = D1
             dinfo_all%VAR_1d(:) = dinfo%VAR_1d(:)

          endselect

       case( 2 )

          D1 = size( dinfo%VAR_2d(:,:), 1 )
          D2 = size( dinfo%VAR_2d(:,:), 2 )

          GD1 = D1 * nprocs_x_out
          GD2 = D2 * nprocs_y_out

          allocate( dinfo_all%VAR_2d( GD1, GD2) )
          dinfo_all%dim_size(1) = GD1
          dinfo_all%dim_size(2) = GD2

          select case( dinfo_all%datatype )
          case( FILE_REAL4 )
             call gather_xy( ismaster,             &
                             SP,                   &
                             nprocs_x_out,         &
                             nprocs_y_out,         &
                             D1, D2,               &
                             dinfo%VAR_2d(:,:),    &
                             dinfo_all%VAR_2d(:,:) )
          case( FILE_REAL8 )
             call gather_xy( ismaster,             &
                             DP,                   &
                             nprocs_x_out,         &
                             nprocs_y_out,         &
                             D1, D2,               &
                             dinfo%VAR_2d(:,:),    &
                             dinfo_all%VAR_2d(:,:) )
          endselect

       case( 3 )

          D1 = size( dinfo%VAR_3d(:,:,:), 1 )
          D2 = size( dinfo%VAR_3d(:,:,:), 2 )
          D3 = size( dinfo%VAR_3d(:,:,:), 3 )

          if ( dinfo_all%transpose ) then

             GD1 = D1
             GD2 = D2 * nprocs_x_out
             GD3 = D3 * nprocs_y_out

             allocate( dinfo_all%VAR_3d( GD1, GD2, GD3 ) )
             dinfo_all%dim_size(1) = GD1
             dinfo_all%dim_size(2) = GD2
             dinfo_all%dim_size(3) = GD3

             select case( dinfo_all%datatype )
             case( FILE_REAL4 )
                call gather_zxy( ismaster,               &
                                 SP,                     &
                                 nprocs_x_out,           &
                                 nprocs_y_out,           &
                                 D1, D2, D3,             &
                                 dinfo%VAR_3d(:,:,:),    &
                                 dinfo_all%VAR_3d(:,:,:) )
             case( FILE_REAL8 )
                call gather_zxy( ismaster,               &
                                 DP,                     &
                                 nprocs_x_out,           &
                                 nprocs_y_out,           &
                                 D1, D2, D3,             &
                                 dinfo%VAR_3d(:,:,:),    &
                                 dinfo_all%VAR_3d(:,:,:) )
             endselect

          else

             GD1 = D1 * nprocs_x_out
             GD2 = D2 * nprocs_y_out
             GD3 = D3

             allocate( dinfo_all%VAR_3d( GD1, GD2, GD3 ) )
             dinfo_all%dim_size(1) = GD1
             dinfo_all%dim_size(2) = GD2
             dinfo_all%dim_size(3) = GD3

             select case( dinfo_all%datatype )
             case( FILE_REAL4 )
                call gather_xyz( ismaster,               &
                                 SP,                     &
                                 nprocs_x_out,           &
                                 nprocs_y_out,           &
                                 D1, D2, D3,             &
                                 dinfo%VAR_3d(:,:,:),    &
                                 dinfo_all%VAR_3d(:,:,:) )
             case( FILE_REAL8 )
                call gather_xyz( ismaster,               &
                                 DP,                     &
                                 nprocs_x_out,           &
                                 nprocs_y_out,           &
                                 D1, D2, D3,             &
                                 dinfo%VAR_3d(:,:,:),    &
                                 dinfo_all%VAR_3d(:,:,:) )
             endselect

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
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1,           &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(RP) :: recv(D1,nprocs_x_out*nprocs_y_out)

    integer  :: i
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHER( send_SP(:),           &
                        D1,                   &
                        MPI_REAL,             &
                        recv(:,:),            &
                        D1,                   &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHER( send_DP(:),           &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:),            &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do i = 1, D1
             dout( i + (px-1) * D1 ) = recv(i,p)
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_x

  !-----------------------------------------------------------------------------
  subroutine gather_u( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1,           &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(RP) :: recv(D1,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: iloc

    integer  :: i
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D1
       else
          recvicnt(px) = D1 - 1
       endif

       recvpcnt(p) = recvicnt(px)
       recvploc(p) = D1 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHERV( send_SP(:),           &
                         D1,                   &
                         MPI_REAL,             &
                         recv(:,:),            &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHERV( send_DP(:),           &
                         D1,                   &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:),            &
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
             iloc = sum( recvicnt(1:px-1) )
          endif

          ! rearrangement of data array
          do i = 1, recvicnt(px)
             dout( i + iloc ) = recv(i,p)
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_u

  !-----------------------------------------------------------------------------
  subroutine gather_y( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1,           &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(RP) :: recv(D1,nprocs_x_out*nprocs_y_out)

    integer  :: j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHER( send_SP(:),           &
                        D1,                   &
                        MPI_REAL,             &
                        recv(:,:),            &
                        D1,                   &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHER( send_DP(:),           &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:),            &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do j = 1, D1
             dout( j + (py-1) * D1 ) = recv(j,p)
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_y

  !-----------------------------------------------------------------------------
  subroutine gather_v( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1,           &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(RP) :: recv(D1,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: jloc

    integer  :: j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( py == 1 ) then
          recvjcnt(py) = D1
       else
          recvjcnt(py) = D1 - 1
       endif

       recvpcnt(p) = recvjcnt(py)
       recvploc(p) = D1 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHERV( send_SP(:),           &
                         D1,                   &
                         MPI_REAL,             &
                         recv(:,:),            &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHERV( send_DP(:),           &
                         D1,                   &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:),            &
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
             jloc = sum( recvjcnt(1:py-1) )
          endif

          ! rearrangement of data array
          do j = 1, recvjcnt(py)
             dout( j + jloc ) = recv(j,p)
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_v

  !-----------------------------------------------------------------------------
  subroutine gather_x_halo( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       xhalo,        &
       D1,           &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: xhalo
    integer,  intent(in)  :: D1

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(RP) :: recv(D1,nprocs_x_out*nprocs_y_out)

    integer  :: i
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHER( send_SP(:),           &
                        D1,                   &
                        MPI_REAL,             &
                        recv(:,:),            &
                        D1,                   &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHER( send_DP(:),           &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:),            &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do i = 1, D1
             dout( i + (px-1) * (D1-2*xhalo) ) = recv(i,p)
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_x_halo

  !-----------------------------------------------------------------------------
  subroutine gather_u_halo( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       xhalo,        &
       D1,           &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: xhalo
    integer,  intent(in)  :: D1

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(RP) :: recv(D1,nprocs_x_out*nprocs_y_out)

    integer  :: i
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHER( send_SP(:),           &
                        D1,                   &
                        MPI_REAL,             &
                        recv(:,:),            &
                        D1,                   &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHER( send_DP(:),           &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:),            &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do i = 1, D1
             dout( i + (px-1) * (D1-2*xhalo-1) ) = recv(i,p)
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_u_halo

  !-----------------------------------------------------------------------------
  subroutine gather_y_halo( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       yhalo,        &
       D1,           &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: yhalo
    integer,  intent(in)  :: D1

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(RP) :: recv(D1,nprocs_x_out*nprocs_y_out)

    integer  :: j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHER( send_SP(:),           &
                        D1,                   &
                        MPI_REAL,             &
                        recv(:,:),            &
                        D1,                   &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHER( send_DP(:),           &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:),            &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do j = 1, D1
             dout( j + (py-1) * (D1-2*yhalo) ) = recv(j,p)
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_y_halo

  !-----------------------------------------------------------------------------
  subroutine gather_v_halo( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       yhalo,        &
       D1,           &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: yhalo
    integer,  intent(in)  :: D1

    real(RP), intent(in)  :: din (:)
    real(RP), intent(out) :: dout(:)

    real(SP) :: send_SP(D1)
    real(DP) :: send_DP(D1)

    real(RP) :: recv(D1,nprocs_x_out*nprocs_y_out)

    integer  :: j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:) = real( din(:), kind=SP )

       call MPI_GATHER( send_SP(:),           &
                        D1,                   &
                        MPI_REAL,             &
                        recv(:,:),            &
                        D1,                   &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:) = real( din(:), kind=DP )

       call MPI_GATHER( send_DP(:),           &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:),            &
                        D1,                   &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do j = 1, D1
             dout( j + (py-1) * (D1-2*yhalo-1) ) = recv(j,p)
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_v_halo

  !-----------------------------------------------------------------------------
  subroutine gather_xy( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2,       &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(RP) :: recv(D1,D2,nprocs_x_out*nprocs_y_out)

    integer  :: i, j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHER( send_SP(:,:),         &
                        D1 * D2,              &
                        MPI_REAL,             &
                        recv(:,:,:),          &
                        D1 * D2,              &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHER( send_DP(:,:),         &
                        D1 * D2,              &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:,:),          &
                        D1 * D2,              &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do j = 1, D2
          do i = 1, D1
             dout( i + (px-1) * D1, &
                   j + (py-1) * D2  ) = recv(i,j,p)
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_xy

  !-----------------------------------------------------------------------------
  subroutine gather_uy( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2,       &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(RP) :: recv(D1,D2,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: iloc

    integer  :: i, j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D1
       else
          recvicnt(px) = D1 - 1
       endif

       recvpcnt(p) = recvicnt(px) * D2
       recvploc(p) = D1 * D2 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:),         &
                         D1 * D2,              &
                         MPI_REAL,             &
                         recv(:,:,:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:),         &
                         D1 * D2,              &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:),          &
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
             iloc = sum( recvicnt(1:px-1) )
          endif

          ! rearrangement of data array
          do j = 1, D2
          do i = 1, recvicnt(px)
             dout( i + iloc,        &
                   j + (py-1) * D2  ) = recv(i,j,p)
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_uy

  !-----------------------------------------------------------------------------
  subroutine gather_xv( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2,       &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(RP) :: recv(D1,D2,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: jloc

    integer  :: i, j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( py == 1 ) then
          recvjcnt(py) = D2
       else
          recvjcnt(py) = D2 - 1
       endif

       recvpcnt(p) = D1 * recvjcnt(py)
       recvploc(p) = D1 * D2 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:),         &
                         D1 * D2,              &
                         MPI_REAL,             &
                         recv(:,:,:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:),         &
                         D1 * D2,              &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:),          &
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
             jloc = sum( recvjcnt(1:py-1) )
          endif

          ! rearrangement of data array
          do j = 1, recvjcnt(py)
          do i = 1, D1
             dout( i + (px-1) * D1, &
                   j + jloc         ) = recv(i,j,p)
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_xv

  !-----------------------------------------------------------------------------
  subroutine gather_uv( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2,       &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(RP) :: recv(D1,D2,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: iloc, jloc

    integer  :: i, j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D1
       else
          recvicnt(px) = D1 - 1
       endif
       if ( py == 1 ) then
          recvjcnt(py) = D2
       else
          recvjcnt(py) = D2 - 1
       endif

       recvpcnt(p) = recvicnt(px) * recvjcnt(py)
       recvploc(p) = D1 * D2 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:),         &
                         D1 * D2,              &
                         MPI_REAL,             &
                         recv(:,:,:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:),         &
                         D1 * D2,              &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:),          &
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
             iloc = sum( recvicnt(1:px-1) )
          endif
          if ( py == 1 ) then
             jloc = 0
          else
             jloc = sum( recvjcnt(1:py-1) )
          endif

          ! rearrangement of data array
          do j = 1, recvjcnt(py)
          do i = 1, recvicnt(px)
             dout( i + iloc, &
                   j + jloc  ) = recv(i,j,p)
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_uv

  !-----------------------------------------------------------------------------
  subroutine gather_nx( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2,       &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(RP) :: recv(D1,D2,nprocs_x_out*nprocs_y_out)

    integer  :: v, i
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHER( send_SP(:,:),         &
                        D1 * D2,              &
                        MPI_REAL,             &
                        recv(:,:,:),          &
                        D1 * D2,              &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHER( send_DP(:,:),         &
                        D1 * D2,              &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:,:),          &
                        D1 * D2,              &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do i = 1, D2
          do v = 1, D1
             dout( v,              &
                   i + (px-1) * D2 ) = recv(v,i,p)
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_nx

  !-----------------------------------------------------------------------------
  subroutine gather_nu( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2,       &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(RP) :: recv(D1,D2,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: iloc

    integer  :: v, i
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D2
       else
          recvicnt(px) = D2 - 1
       endif

       recvpcnt(p) = D1 * recvicnt(px)
       recvploc(p) = D1 * D2 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:),         &
                         D1 * D2,              &
                         MPI_REAL,             &
                         recv(:,:,:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:),         &
                         D1 * D2,              &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:),          &
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
             iloc = sum( recvicnt(1:px-1) )
          endif

          ! rearrangement of data array
          do i = 1, recvicnt(px)
          do v = 1, D1
             dout( v,       &
                   i + iloc ) = recv(v,i,p)
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_nu

  !-----------------------------------------------------------------------------
  subroutine gather_ny( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2,       &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(RP) :: recv(D1,D2,nprocs_x_out*nprocs_y_out)

    integer  :: v, j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHER( send_SP(:,:),         &
                        D1 * D2,              &
                        MPI_REAL,             &
                        recv(:,:,:),          &
                        D1 * D2,              &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHER( send_DP(:,:),         &
                        D1 * D2,              &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:,:),          &
                        D1 * D2,              &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do j = 1, D2
          do v = 1, D1
             dout( v,              &
                   j + (py-1) * D2 ) = recv(v,j,p)
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_ny

  !-----------------------------------------------------------------------------
  subroutine gather_nv( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2,       &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2

    real(RP), intent(in)  :: din (:,:)
    real(RP), intent(out) :: dout(:,:)

    real(SP) :: send_SP(D1,D2)
    real(DP) :: send_DP(D1,D2)

    real(RP) :: recv(D1,D2,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: jloc

    integer  :: v, j
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( py == 1 ) then
          recvjcnt(py) = D2
       else
          recvjcnt(py) = D2 - 1
       endif

       recvpcnt(p) = D1 * recvjcnt(py)
       recvploc(p) = D1 * D2 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:) = real( din(:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:),         &
                         D1 * D2,              &
                         MPI_REAL,             &
                         recv(:,:,:),          &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:) = real( din(:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:),         &
                         D1 * D2,              &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:),          &
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
             jloc = sum( recvjcnt(1:py-1) )
          endif

          ! rearrangement of data array
          do j = 1, recvjcnt(py)
          do v = 1, D1
             dout( v,       &
                   j + jloc ) = recv(v,j,p)
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_nv

  !-----------------------------------------------------------------------------
  subroutine gather_zxy( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2, D3,   &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2
    integer,  intent(in)  :: D3

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(RP) :: recv(D1,D2,D3,nprocs_x_out*nprocs_y_out)

    integer  :: i, j, k
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHER( send_SP(:,:,:),       &
                        D1 * D2 * D3,         &
                        MPI_REAL,             &
                        recv(:,:,:,:),        &
                        D1 * D2 * D3,         &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHER( send_DP(:,:,:),       &
                        D1 * D2 * D3,         &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:,:,:),        &
                        D1 * D2 * D3,         &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do j = 1, D3
          do i = 1, D2
          do k = 1, D1
             dout( k,               &
                   i + (px-1) * D2, &
                   j + (py-1) * D3  ) = recv(k,i,j,p)
          enddo
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_zxy

  !-----------------------------------------------------------------------------
  subroutine gather_zuy( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2, D3,   &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2
    integer,  intent(in)  :: D3

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(RP) :: recv(D1,D2,D3,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: iloc

    integer  :: i, j, k
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D2
       else
          recvicnt(px) = D2 - 1
       endif

       recvpcnt(p) = D1 * recvicnt(px) * D3
       recvploc(p) = D1 * D2 * D3 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_REAL,             &
                         recv(:,:,:,:),        &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:,:),        &
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
             iloc = sum( recvicnt(1:px-1) )
          endif

          ! rearrangement of data array
          do j = 1, D3
          do i = 1, recvicnt(px)
          do k = 1, D1
             dout( k,               &
                   i + iloc,        &
                   j + (py-1) * D3  ) = recv(k,i,j,p)
          enddo
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_zuy

  !-----------------------------------------------------------------------------
  subroutine gather_zxv( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2, D3,   &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2
    integer,  intent(in)  :: D3

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(RP) :: recv(D1,D2,D3,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: jloc

    integer  :: i, j, k
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( py == 1 ) then
          recvjcnt(py) = D3
       else
          recvjcnt(py) = D3 - 1
       endif

       recvpcnt(p) = D1 * D2 * recvjcnt(py)
       recvploc(p) = D1 * D2 * D3 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_REAL,             &
                         recv(:,:,:,:),        &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:,:),        &
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
             jloc = sum( recvjcnt(1:py-1) )
          endif

          ! rearrangement of data array
          do j = 1, recvjcnt(py)
          do i = 1, D2
          do k = 1, D1
             dout( k,               &
                   i + (px-1) * D2, &
                   j + jloc         ) = recv(k,i,j,p)
          enddo
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_zxv

  !-----------------------------------------------------------------------------
  subroutine gather_zuv( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2, D3,   &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2
    integer,  intent(in)  :: D3

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(RP) :: recv(D1,D2,D3,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: iloc, jloc

    integer  :: i, j, k
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D2
       else
          recvicnt(px) = D2 - 1
       endif
       if ( py == 1 ) then
          recvjcnt(py) = D3
       else
          recvjcnt(py) = D3 - 1
       endif

       recvpcnt(p) = D1 * recvicnt(px) * recvjcnt(py)
       recvploc(p) = D1 * D2 * D3 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_REAL,             &
                         recv(:,:,:,:),        &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:,:),        &
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
             iloc = sum( recvicnt(1:px-1) )
          endif
          if ( py == 1 ) then
             jloc = 0
          else
             jloc = sum( recvjcnt(1:py-1) )
          endif

          ! rearrangement of data array
          do j = 1, recvjcnt(py)
          do i = 1, recvicnt(px)
          do k = 1, D1
             dout( k,        &
                   i + iloc, &
                   j + jloc  ) = recv(k,i,j,p)
          enddo
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_zuv

  !-----------------------------------------------------------------------------
  subroutine gather_xyz( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2, D3,   &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2
    integer,  intent(in)  :: D3

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(RP) :: recv(D1,D2,D3,nprocs_x_out*nprocs_y_out)

    integer  :: i, j, k
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    select case( MPI_RP )
    case( SP )
       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHER( send_SP(:,:,:),       &
                        D1 * D2 * D3,         &
                        MPI_REAL,             &
                        recv(:,:,:,:),        &
                        D1 * D2 * D3,         &
                        MPI_REAL,             &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    case( DP )
       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHER( send_DP(:,:,:),       &
                        D1 * D2 * D3,         &
                        MPI_DOUBLE_PRECISION, &
                        recv(:,:,:,:),        &
                        D1 * D2 * D3,         &
                        MPI_DOUBLE_PRECISION, &
                        PRC_masterrank,       &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    endselect

    if ( ismaster ) then

       p = 1
       do py = 1, nprocs_y_out
       do px = 1, nprocs_x_out
          ! rearrangement of data array
          do k = 1, D3
          do j = 1, D2
          do i = 1, D1
             dout( i + (px-1) * D1, &
                   j + (py-1) * D2, &
                   k                ) = recv(i,j,k,p)
          enddo
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_xyz

  !-----------------------------------------------------------------------------
  subroutine gather_uyz( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2, D3,   &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2
    integer,  intent(in)  :: D3

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(RP) :: recv(D1,D2,D3,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: iloc

    integer  :: i, j, k
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D1
       else
          recvicnt(px) = D1 - 1
       endif

       recvpcnt(p) = recvicnt(px) * D2 * D3
       recvploc(p) = D1 * D2 * D3 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_REAL,             &
                         recv(:,:,:,:),        &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:,:),        &
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
             iloc = sum( recvicnt(1:px-1) )
          endif

          ! rearrangement of data array
          do k = 1, D3
          do j = 1, D2
          do i = 1, recvicnt(px)
             dout( i + iloc,        &
                   j + (py-1) * D2, &
                   k                ) = recv(i,j,k,p)

          enddo
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_uyz

  !-----------------------------------------------------------------------------
  subroutine gather_xvz( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2, D3,   &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2
    integer,  intent(in)  :: D3

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(RP) :: recv(D1,D2,D3,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: jloc

    integer  :: i, j, k
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( py == 1 ) then
          recvjcnt(py) = D2
       else
          recvjcnt(py) = D2 - 1
       endif

       recvpcnt(p) = D1 * recvjcnt(py) * D3
       recvploc(p) = D1 * D2 * D3 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_REAL,             &
                         recv(:,:,:,:),        &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:,:),        &
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
             jloc = sum( recvjcnt(1:py-1) )
          endif

          ! rearrangement of data array
          do k = 1, D3
          do j = 1, recvjcnt(py)
          do i = 1, D1
             dout( i + (px-1) * D1, &
                   j + jloc,        &
                   k                ) = recv(i,j,k,p)
          enddo
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_xvz

  !-----------------------------------------------------------------------------
  subroutine gather_uvz( &
       ismaster,     &
       MPI_RP,       &
       nprocs_x_out, &
       nprocs_y_out, &
       D1, D2, D3,   &
       din,          &
       dout          )
    use mpi
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    implicit none

    logical,  intent(in)  :: ismaster
    integer,  intent(in)  :: MPI_RP
    integer,  intent(in)  :: nprocs_x_out
    integer,  intent(in)  :: nprocs_y_out
    integer,  intent(in)  :: D1
    integer,  intent(in)  :: D2
    integer,  intent(in)  :: D3

    real(RP), intent(in)  :: din (:,:,:)
    real(RP), intent(out) :: dout(:,:,:)

    real(SP) :: send_SP(D1,D2,D3)
    real(DP) :: send_DP(D1,D2,D3)

    real(RP) :: recv(D1,D2,D3,nprocs_x_out*nprocs_y_out)

    integer  :: recvpcnt(nprocs_x_out*nprocs_y_out)
    integer  :: recvploc(nprocs_x_out*nprocs_y_out)
    integer  :: recvicnt(nprocs_x_out)
    integer  :: recvjcnt(nprocs_y_out)
    integer  :: iloc, jloc

    integer  :: i, j, k
    integer  :: p, px, py
    integer  :: ierr
    !---------------------------------------------------------------------------

    p = 1
    do py = 1, nprocs_y_out
    do px = 1, nprocs_x_out
       if ( px == 1 ) then
          recvicnt(px) = D1
       else
          recvicnt(px) = D1 - 1
       endif
       if ( py == 1 ) then
          recvjcnt(py) = D2
       else
          recvjcnt(py) = D2 - 1
       endif

       recvpcnt(p) = recvicnt(px) * recvjcnt(py) * D3
       recvploc(p) = D1 * D2 * D3 * (p-1)

       p = p + 1
    enddo
    enddo

    select case( MPI_RP )
    case( SP )
       send_SP(:,:,:) = real( din(:,:,:), kind=SP )

       call MPI_GATHERV( send_SP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_REAL,             &
                         recv(:,:,:,:),        &
                         recvpcnt(:),          &
                         recvploc(:),          &
                         MPI_REAL,             &
                         PRC_masterrank,       &
                         PRC_LOCAL_COMM_WORLD, &
                         ierr                  )
    case( DP )
       send_DP(:,:,:) = real( din(:,:,:), kind=DP )

       call MPI_GATHERV( send_DP(:,:,:),       &
                         D1 * D2 * D3,         &
                         MPI_DOUBLE_PRECISION, &
                         recv(:,:,:,:),        &
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
             iloc = sum( recvicnt(1:px-1) )
          endif
          if ( py == 1 ) then
             jloc = 0
          else
             jloc = sum( recvjcnt(1:py-1) )
          endif

          ! rearrangement of data array
          do k = 1, D3
          do j = 1, recvjcnt(py)
          do i = 1, recvicnt(px)
             dout( i + iloc, &
                   j + jloc, &
                   k         ) = recv(i,j,k,p)
          enddo
          enddo
          enddo
          p = p + 1
       enddo
       enddo

    endif

    return
  end subroutine gather_uvz

end module mod_sno_comm
