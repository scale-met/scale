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
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNO_comm_globalaxis( &
       output_single, &
       nprocs_x_out,  &
       nprocs_y_out,  &
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
       axisinfo
    implicit none

    logical,          intent(in)    :: output_single                         ! output single file when using MPI?
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (input)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (input)
    integer,          intent(in)    :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(in)    :: ainfo(:)                              ! axis information                   (input)
    type(axisinfo),   intent(out)   :: ainfo_all(:)                          ! axis information                   (output)

    real(DP), allocatable :: send_1d_DP(:)
    real(DP), allocatable :: send_2d_DP(:,:)
    real(DP), allocatable :: send_3d_DP(:,:,:)

    real(DP), allocatable :: recv_1d_DP(:,:)
    real(DP), allocatable :: recv_2d_DP(:,:,:)
    real(DP), allocatable :: recv_3d_DP(:,:,:,:)

    real(SP), allocatable :: send_1d_SP(:)
    real(SP), allocatable :: send_2d_SP(:,:)
    real(SP), allocatable :: send_3d_SP(:,:,:)

    real(SP), allocatable :: recv_1d_SP(:,:)
    real(SP), allocatable :: recv_2d_SP(:,:,:)
    real(SP), allocatable :: recv_3d_SP(:,:,:,:)

    integer  :: i, j, k, n, v, t
    integer  :: p, px, py
    integer  :: D1, D2, D3
    integer  :: nprocs

    integer  :: sendcount, recvcount
    integer  :: datatype
    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( output_single ) then

       nprocs = nprocs_x_out * nprocs_y_out

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

             sendcount = D1
             recvcount = D1

             select case( ainfo_all(n)%varname )
             case('x','xh','CX','FX','CDX','FDX','CBFX','FBFX')

                allocate( ainfo_all(n)%AXIS_1d( D1*nprocs_x_out ) )

                ainfo_all(n)%dim_size(1) = D1 * nprocs_x_out

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )

                   datatype = MPI_REAL

                   allocate( send_1d_SP(D1       ) )
                   allocate( recv_1d_SP(D1,nprocs) )

                   send_1d_SP(:) = real( ainfo(n)%AXIS_1d(:), kind=SP )

                   call MPI_GATHER( send_1d_SP(:),        &
                                    sendcount,            &
                                    datatype,             &
                                    recv_1d_SP(:,:),      &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do i = 1, D1
                         ainfo_all(n)%AXIS_1d( i + (px-1) * D1 ) = recv_1d_SP(i,p)
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_1d_SP )
                   deallocate( recv_1d_SP )

                case( FILE_REAL8 )

                   datatype = MPI_DOUBLE_PRECISION

                   allocate( send_1d_DP(D1       ) )
                   allocate( recv_1d_DP(D1,nprocs) )

                   send_1d_DP(:) = real( ainfo(n)%AXIS_1d(:), kind=DP )

                   call MPI_GATHER( send_1d_DP(:),        &
                                    sendcount,            &
                                    datatype,             &
                                    recv_1d_DP(:,:),      &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do i = 1, D1
                         ainfo_all(n)%AXIS_1d( i + (px-1) * D1 ) = recv_1d_DP(i,p)
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_1d_DP )
                   deallocate( recv_1d_DP )

                endselect

             case('y','yh','CY','FY','CDY','FDY','CBFY','FBFY')

                allocate( ainfo_all(n)%AXIS_1d( D1*nprocs_y_out ) )

                ainfo_all(n)%dim_size(1) = D1 * nprocs_y_out

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )

                   datatype = MPI_REAL

                   allocate( send_1d_SP(D1       ) )
                   allocate( recv_1d_SP(D1,nprocs) )

                   send_1d_SP(:) = real( ainfo(n)%AXIS_1d(:), kind=SP )

                   call MPI_GATHER( send_1d_SP(:),        &
                                    sendcount,            &
                                    datatype,             &
                                    recv_1d_SP(:,:),      &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do j = 1, D1
                         ainfo_all(n)%AXIS_1d( j + (py-1) * D1 ) = recv_1d_SP(j,p)
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_1d_SP )
                   deallocate( recv_1d_SP )

                case( FILE_REAL8 )

                   datatype = MPI_DOUBLE_PRECISION

                   allocate( send_1d_DP(D1       ) )
                   allocate( recv_1d_DP(D1,nprocs) )

                   send_1d_DP(:) = real( ainfo(n)%AXIS_1d(:), kind=DP )

                   call MPI_GATHER( send_1d_DP(:),        &
                                    sendcount,            &
                                    datatype,             &
                                    recv_1d_DP(:,:),      &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do j = 1, D1
                         ainfo_all(n)%AXIS_1d( j + (py-1) * D1 ) = recv_1d_DP(j,p)
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_1d_DP )
                   deallocate( recv_1d_DP )

                endselect

             case default

                allocate( ainfo_all(n)%AXIS_1d( D1 ) )

                ainfo_all(n)%dim_size(1) = D1
                ainfo_all(n)%AXIS_1d(:) = ainfo(n)%AXIS_1d(:)

             endselect

          case( 2 )

             D1 = size( ainfo(n)%AXIS_2d(:,:), 1 )
             D2 = size( ainfo(n)%AXIS_2d(:,:), 2 )

             sendcount = D1 * D2
             recvcount = D1 * D2

             select case( index( ainfo_all(n)%varname, '_bnds' ) )
             case( 0 )

                allocate( ainfo_all(n)%AXIS_2d( D1*nprocs_x_out, D2*nprocs_y_out ) )

                ainfo_all(n)%dim_size(1) = D1 * nprocs_x_out
                ainfo_all(n)%dim_size(2) = D2 * nprocs_y_out

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )

                   datatype = MPI_REAL

                   allocate( send_2d_SP(D1,D2       ) )
                   allocate( recv_2d_SP(D1,D2,nprocs) )

                   send_2d_SP(:,:) = real( ainfo(n)%AXIS_2d(:,:), kind=SP )

                   call MPI_GATHER( send_2d_SP(:,:),      &
                                    sendcount,            &
                                    datatype,             &
                                    recv_2d_SP(:,:,:),    &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do j = 1, D2
                      do i = 1, D1
                         ainfo_all(n)%AXIS_2d( i + (px-1) * D1, &
                                               j + (py-1) * D2  ) = recv_2d_SP(i,j,p)
                      enddo
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_2d_SP )
                   deallocate( recv_2d_SP )

                case( FILE_REAL8 )

                   datatype = MPI_DOUBLE_PRECISION

                   allocate( send_2d_DP(D1,D2       ) )
                   allocate( recv_2d_DP(D1,D2,nprocs) )

                   send_2d_DP(:,:) = real( ainfo(n)%AXIS_2d(:,:), kind=DP )

                   call MPI_GATHER( send_2d_DP(:,:),      &
                                    sendcount,            &
                                    datatype,             &
                                    recv_2d_DP(:,:,:),    &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do j = 1, D2
                      do i = 1, D1
                         ainfo_all(n)%AXIS_2d( i + (px-1) * D1, &
                                               j + (py-1) * D2  ) = recv_2d_DP(i,j,p)
                      enddo
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_2d_DP )
                   deallocate( recv_2d_DP )

                endselect

             case default

                select case( ainfo_all(n)%varname )
                case('x_bnds','xh_bnds')

                   allocate( ainfo_all(n)%AXIS_2d( D1, D2*nprocs_x_out ) )

                   ainfo_all(n)%dim_size(1) = D1
                   ainfo_all(n)%dim_size(2) = D2 * nprocs_x_out

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )

                      datatype = MPI_REAL

                      allocate( send_2d_SP(D1,D2       ) )
                      allocate( recv_2d_SP(D1,D2,nprocs) )

                      send_2d_SP(:,:) = real( ainfo(n)%AXIS_2d(:,:), kind=SP )

                      call MPI_GATHER( send_2d_SP(:,:),      &
                                       sendcount,            &
                                       datatype,             &
                                       recv_2d_SP(:,:,:),    &
                                       recvcount,            &
                                       datatype,             &
                                       PRC_masterrank,       &
                                       PRC_LOCAL_COMM_WORLD, &
                                       ierr                  )

                      p = 1
                      do py = 1, nprocs_y_out
                      do px = 1, nprocs_x_out
                         ! rearrangement of data array
                         do i = 1, D2
                         do v = 1, D1
                            ainfo_all(n)%AXIS_2d( v,              &
                                                  i + (px-1) * D2 ) = recv_2d_SP(v,i,p)
                         enddo
                         enddo
                         p = p + 1
                      enddo
                      enddo

                      deallocate( send_2d_SP )
                      deallocate( recv_2d_SP )

                   case( FILE_REAL8 )

                      datatype = MPI_DOUBLE_PRECISION

                      allocate( send_2d_DP(D1,D2       ) )
                      allocate( recv_2d_DP(D1,D2,nprocs) )

                      send_2d_DP(:,:) = real( ainfo(n)%AXIS_2d(:,:), kind=DP )

                      call MPI_GATHER( send_2d_DP(:,:),      &
                                       sendcount,            &
                                       datatype,             &
                                       recv_2d_DP(:,:,:),    &
                                       recvcount,            &
                                       datatype,             &
                                       PRC_masterrank,       &
                                       PRC_LOCAL_COMM_WORLD, &
                                       ierr                  )

                      p = 1
                      do py = 1, nprocs_y_out
                      do px = 1, nprocs_x_out
                         ! rearrangement of data array
                         do i = 1, D2
                         do v = 1, D1
                            ainfo_all(n)%AXIS_2d( v,              &
                                                  i + (px-1) * D2 ) = recv_2d_DP(v,i,p)
                         enddo
                         enddo
                         p = p + 1
                      enddo
                      enddo

                      deallocate( send_2d_DP )
                      deallocate( recv_2d_DP )

                   endselect

                case('y_bnds','yh_bnds')

                   allocate( ainfo_all(n)%AXIS_2d( D1, D2*nprocs_y_out ) )

                   ainfo_all(n)%dim_size(1) = D1
                   ainfo_all(n)%dim_size(2) = D2 * nprocs_y_out

                   select case( ainfo_all(n)%datatype )
                   case( FILE_REAL4 )

                      datatype = MPI_REAL

                      allocate( send_2d_SP(D1,D2       ) )
                      allocate( recv_2d_SP(D1,D2,nprocs) )

                      send_2d_SP(:,:) = real( ainfo(n)%AXIS_2d(:,:), kind=SP )

                      call MPI_GATHER( send_2d_SP(:,:),      &
                                       sendcount,            &
                                       datatype,             &
                                       recv_2d_SP(:,:,:),    &
                                       recvcount,            &
                                       datatype,             &
                                       PRC_masterrank,       &
                                       PRC_LOCAL_COMM_WORLD, &
                                       ierr                  )

                      p = 1
                      do py = 1, nprocs_y_out
                      do px = 1, nprocs_x_out
                         ! rearrangement of data array
                         do j = 1, D2
                         do v = 1, D1
                            ainfo_all(n)%AXIS_2d( v,              &
                                                  j + (py-1) * D2 ) = recv_2d_SP(v,j,p)
                         enddo
                         enddo
                         p = p + 1
                      enddo
                      enddo

                      deallocate( send_2d_SP )
                      deallocate( recv_2d_SP )

                   case( FILE_REAL8 )

                      datatype = MPI_DOUBLE_PRECISION

                      allocate( send_2d_DP(D1,D2       ) )
                      allocate( recv_2d_DP(D1,D2,nprocs) )

                      send_2d_DP(:,:) = real( ainfo(n)%AXIS_2d(:,:), kind=DP )

                      call MPI_GATHER( send_2d_DP(:,:),      &
                                       sendcount,            &
                                       datatype,             &
                                       recv_2d_DP(:,:,:),    &
                                       recvcount,            &
                                       datatype,             &
                                       PRC_masterrank,       &
                                       PRC_LOCAL_COMM_WORLD, &
                                       ierr                  )

                      p = 1
                      do py = 1, nprocs_y_out
                      do px = 1, nprocs_x_out
                         ! rearrangement of data array
                         do j = 1, D2
                         do v = 1, D1
                            ainfo_all(n)%AXIS_2d( v,              &
                                                  j + (py-1) * D2 ) = recv_2d_DP(v,j,p)
                         enddo
                         enddo
                         p = p + 1
                      enddo
                      enddo

                      deallocate( send_2d_DP )
                      deallocate( recv_2d_DP )

                   endselect

                case default

                   allocate( ainfo_all(n)%AXIS_2d( D1, D2 ) )

                   ainfo_all(n)%dim_size(1) = D1
                   ainfo_all(n)%dim_size(2) = D2
                   ainfo_all(n)%AXIS_2d(:,:) = ainfo(n)%AXIS_2d(:,:)

                endselect

             endselect

          case( 3 )

             D1 = size( ainfo(n)%AXIS_3d(:,:,:), 1 )
             D2 = size( ainfo(n)%AXIS_3d(:,:,:), 2 )
             D3 = size( ainfo(n)%AXIS_3d(:,:,:), 3 )

             sendcount = D1 * D2 * D3
             recvcount = D1 * D2 * D3

             if ( ainfo_all(n)%transpose ) then

                allocate( ainfo_all(n)%AXIS_3d( D1, D2*nprocs_x_out, D3*nprocs_y_out ) )

                ainfo_all(n)%dim_size(1) = D1
                ainfo_all(n)%dim_size(2) = D2 * nprocs_x_out
                ainfo_all(n)%dim_size(3) = D3 * nprocs_y_out

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )

                   datatype = MPI_REAL

                   allocate( send_3d_SP(D1,D2,D3       ) )
                   allocate( recv_3d_SP(D1,D2,D3,nprocs) )

                   send_3d_SP(:,:,:) = real( ainfo(n)%AXIS_3d(:,:,:), kind=SP )

                   call MPI_GATHER( send_3d_SP(:,:,:),    &
                                    sendcount,            &
                                    datatype,             &
                                    recv_3d_SP(:,:,:,:),  &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do j = 1, D3
                      do i = 1, D2
                      do k = 1, D1
                         ainfo_all(n)%AXIS_3d( k,               &
                                               i + (px-1) * D2, &
                                               j + (py-1) * D3  ) = recv_3d_SP(k,i,j,p)
                      enddo
                      enddo
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_3d_SP )
                   deallocate( recv_3d_SP )

                case( FILE_REAL8 )

                   datatype = MPI_DOUBLE_PRECISION

                   allocate( send_3d_DP(D1,D2,D3       ) )
                   allocate( recv_3d_DP(D1,D2,D3,nprocs) )

                   send_3d_DP(:,:,:) = real( ainfo(n)%AXIS_3d(:,:,:), kind=DP )

                   call MPI_GATHER( send_3d_DP(:,:,:),    &
                                    sendcount,            &
                                    datatype,             &
                                    recv_3d_DP(:,:,:,:),  &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do j = 1, D3
                      do i = 1, D2
                      do k = 1, D1
                         ainfo_all(n)%AXIS_3d( k,               &
                                               i + (px-1) * D2, &
                                               j + (py-1) * D3  ) = recv_3d_DP(k,i,j,p)
                      enddo
                      enddo
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_3d_DP )
                   deallocate( recv_3d_DP )

                endselect

             else

                allocate( ainfo_all(n)%AXIS_3d( D1*nprocs_x_out, D2*nprocs_y_out, D3 ) )

                ainfo_all(n)%dim_size(1) = D1 * nprocs_x_out
                ainfo_all(n)%dim_size(2) = D2 * nprocs_y_out
                ainfo_all(n)%dim_size(3) = D3

                select case( ainfo_all(n)%datatype )
                case( FILE_REAL4 )

                   datatype = MPI_REAL

                   allocate( send_3d_SP(D1,D2,D3       ) )
                   allocate( recv_3d_SP(D1,D2,D3,nprocs) )

                   send_3d_SP(:,:,:) = real( ainfo(n)%AXIS_3d(:,:,:), kind=SP )

                   call MPI_GATHER( send_3d_SP(:,:,:),    &
                                    sendcount,            &
                                    datatype,             &
                                    recv_3d_SP(:,:,:,:),  &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do k = 1, D3
                      do j = 1, D2
                      do i = 1, D1
                         ainfo_all(n)%AXIS_3d( i + (px-1) * D1, &
                                               j + (py-1) * D2, &
                                               k                ) = recv_3d_SP(i,j,k,p)
                      enddo
                      enddo
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_3d_SP )
                   deallocate( recv_3d_SP )

                case( FILE_REAL8 )

                   datatype = MPI_DOUBLE_PRECISION

                   allocate( send_3d_DP(D1,D2,D3       ) )
                   allocate( recv_3d_DP(D1,D2,D3,nprocs) )

                   send_3d_DP(:,:,:) = real( ainfo(n)%AXIS_3d(:,:,:), kind=DP )

                   call MPI_GATHER( send_3d_DP(:,:,:),    &
                                    sendcount,            &
                                    datatype,             &
                                    recv_3d_DP(:,:,:,:),  &
                                    recvcount,            &
                                    datatype,             &
                                    PRC_masterrank,       &
                                    PRC_LOCAL_COMM_WORLD, &
                                    ierr                  )

                   p = 1
                   do py = 1, nprocs_y_out
                   do px = 1, nprocs_x_out
                      ! rearrangement of data array
                      do k = 1, D3
                      do j = 1, D2
                      do i = 1, D1
                         ainfo_all(n)%AXIS_3d( i + (px-1) * D1, &
                                               j + (py-1) * D2, &
                                               k                ) = recv_3d_DP(i,j,k,p)
                      enddo
                      enddo
                      enddo
                      p = p + 1
                   enddo
                   enddo

                   deallocate( send_3d_DP )
                   deallocate( recv_3d_DP )

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
       output_single, &
       nprocs_x_out,  &
       nprocs_y_out,  &
       dinfo,         &
       dinfo_all      )
    use mpi
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD
    use mod_sno_h, only: &
       iteminfo
    implicit none

    logical,          intent(in)    :: output_single                         ! output single file when using MPI?
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (input)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (input)
    type(iteminfo),   intent(in)    :: dinfo                                 ! variable information               (input)
    type(iteminfo),   intent(out)   :: dinfo_all                             ! variable information               (output)

    real(DP), allocatable :: send_1d_DP(:)
    real(DP), allocatable :: send_2d_DP(:,:)
    real(DP), allocatable :: send_3d_DP(:,:,:)

    real(DP), allocatable :: recv_1d_DP(:,:)
    real(DP), allocatable :: recv_2d_DP(:,:,:)
    real(DP), allocatable :: recv_3d_DP(:,:,:,:)

    real(SP), allocatable :: send_1d_SP(:)
    real(SP), allocatable :: send_2d_SP(:,:)
    real(SP), allocatable :: send_3d_SP(:,:,:)

    real(SP), allocatable :: recv_1d_SP(:,:)
    real(SP), allocatable :: recv_2d_SP(:,:,:)
    real(SP), allocatable :: recv_3d_SP(:,:,:,:)

    integer  :: i, j, k, n, t
    integer  :: p, px, py
    integer  :: D1, D2, D3
    integer  :: nprocs

    integer  :: sendcount, recvcount
    integer  :: datatype
    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( output_single ) then

       nprocs = nprocs_x_out * nprocs_y_out

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

          sendcount = D1
          recvcount = D1

          select case( dinfo_all%dim_name(1) )
          case('x','xh','CX','FX','CDX','FDX','CBFX','FBFX')

             allocate( dinfo_all%VAR_1d( D1*nprocs_x_out ) )

             dinfo_all%dim_size(1) = D1 * nprocs_x_out

             select case( dinfo_all%datatype )
             case( FILE_REAL4 )

                datatype = MPI_REAL

                allocate( send_1d_SP(D1       ) )
                allocate( recv_1d_SP(D1,nprocs) )

                send_1d_SP(:) = real( dinfo%VAR_1d(:), kind=SP )

                call MPI_GATHER( send_1d_SP(:),        &
                                 sendcount,            &
                                 datatype,             &
                                 recv_1d_SP(:,:),      &
                                 recvcount,            &
                                 datatype,             &
                                 PRC_masterrank,       &
                                 PRC_LOCAL_COMM_WORLD, &
                                 ierr                  )

                p = 1
                do py = 1, nprocs_y_out
                do px = 1, nprocs_x_out
                   ! rearrangement of data array
                   do i = 1, D1
                      dinfo_all%VAR_1d( i + (px-1) * D1 ) = recv_1d_SP(i,p)
                   enddo
                   p = p + 1
                enddo
                enddo

                deallocate( send_1d_SP )
                deallocate( recv_1d_SP )

             case( FILE_REAL8 )

                datatype = MPI_DOUBLE_PRECISION

                allocate( send_1d_DP(D1       ) )
                allocate( recv_1d_DP(D1,nprocs) )

                send_1d_DP(:) = real( dinfo%VAR_1d(:), kind=DP )

                call MPI_GATHER( send_1d_DP(:),        &
                                 sendcount,            &
                                 datatype,             &
                                 recv_1d_DP(:,:),      &
                                 recvcount,            &
                                 datatype,             &
                                 PRC_masterrank,       &
                                 PRC_LOCAL_COMM_WORLD, &
                                 ierr                  )

                p = 1
                do py = 1, nprocs_y_out
                do px = 1, nprocs_x_out
                   ! rearrangement of data array
                   do i = 1, D1
                      dinfo_all%VAR_1d( i + (px-1) * D1 ) = recv_1d_DP(i,p)
                   enddo
                   p = p + 1
                enddo
                enddo

                deallocate( send_1d_DP )
                deallocate( recv_1d_DP )

             endselect

          case('y','yh','CY','FY','CDY','FDY','CBFY','FBFY')

             allocate( dinfo_all%VAR_1d( D1*nprocs_y_out ) )

             dinfo_all%dim_size(1) = D1 * nprocs_y_out

             select case( dinfo_all%datatype )
             case( FILE_REAL4 )

                datatype = MPI_REAL

                allocate( send_1d_SP(D1       ) )
                allocate( recv_1d_SP(D1,nprocs) )

                send_1d_SP(:) = real( dinfo%VAR_1d(:), kind=SP )

                call MPI_GATHER( send_1d_SP(:),        &
                                 sendcount,            &
                                 datatype,             &
                                 recv_1d_SP(:,:),      &
                                 recvcount,            &
                                 datatype,             &
                                 PRC_masterrank,       &
                                 PRC_LOCAL_COMM_WORLD, &
                                 ierr                  )

                p = 1
                do py = 1, nprocs_y_out
                do px = 1, nprocs_x_out
                   ! rearrangement of data array
                   do j = 1, D1
                      dinfo_all%VAR_1d( j + (py-1) * D1 ) = recv_1d_SP(j,p)
                   enddo
                   p = p + 1
                enddo
                enddo

                deallocate( send_1d_SP )
                deallocate( recv_1d_SP )

             case( FILE_REAL8 )

                datatype = MPI_DOUBLE_PRECISION

                allocate( send_1d_DP(D1       ) )
                allocate( recv_1d_DP(D1,nprocs) )

                send_1d_DP(:) = real( dinfo%VAR_1d(:), kind=DP )

                call MPI_GATHER( send_1d_DP(:),        &
                                 sendcount,            &
                                 datatype,             &
                                 recv_1d_DP(:,:),      &
                                 recvcount,            &
                                 datatype,             &
                                 PRC_masterrank,       &
                                 PRC_LOCAL_COMM_WORLD, &
                                 ierr                  )

                p = 1
                do py = 1, nprocs_y_out
                do px = 1, nprocs_x_out
                   ! rearrangement of data array
                   do j = 1, D1
                      dinfo_all%VAR_1d( j + (py-1) * D1 ) = recv_1d_DP(j,p)
                   enddo
                   p = p + 1
                enddo
                enddo

                deallocate( send_1d_DP )
                deallocate( recv_1d_DP )

             endselect

          case default

             allocate( dinfo_all%VAR_1d( D1 ) )

             dinfo_all%dim_size(1) = D1
             dinfo_all%VAR_1d(:) = dinfo%VAR_1d(:)

          endselect

       case( 2 )

          D1 = size( dinfo%VAR_2d(:,:), 1 )
          D2 = size( dinfo%VAR_2d(:,:), 2 )

          sendcount = D1 * D2
          recvcount = D1 * D2

          allocate( dinfo_all%VAR_2d( D1*nprocs_x_out, D2*nprocs_y_out ) )

          dinfo_all%dim_size(1) = D1 * nprocs_x_out
          dinfo_all%dim_size(2) = D2 * nprocs_y_out

          select case( dinfo_all%datatype )
          case( FILE_REAL4 )

             datatype = MPI_REAL

             allocate( send_2d_SP(D1,D2       ) )
             allocate( recv_2d_SP(D1,D2,nprocs) )

             send_2d_SP(:,:) = real( dinfo%VAR_2d(:,:), kind=SP )

             call MPI_GATHER( send_2d_SP(:,:),      &
                              sendcount,            &
                              datatype,             &
                              recv_2d_SP(:,:,:),    &
                              recvcount,            &
                              datatype,             &
                              PRC_masterrank,       &
                              PRC_LOCAL_COMM_WORLD, &
                              ierr                  )

             p = 1
             do py = 1, nprocs_y_out
             do px = 1, nprocs_x_out
                ! rearrangement of data array
                do j = 1, D2
                do i = 1, D1
                   dinfo_all%VAR_2d( i + (px-1) * D1, &
                                     j + (py-1) * D2  ) = recv_2d_SP(i,j,p)
                enddo
                enddo
                p = p + 1
             enddo
             enddo

             deallocate( send_2d_SP )
             deallocate( recv_2d_SP )

          case( FILE_REAL8 )

             datatype = MPI_DOUBLE_PRECISION

             allocate( send_2d_DP(D1,D2       ) )
             allocate( recv_2d_DP(D1,D2,nprocs) )

             send_2d_DP(:,:) = real( dinfo%VAR_2d(:,:), kind=DP )

             call MPI_GATHER( send_2d_DP(:,:),      &
                              sendcount,            &
                              datatype,             &
                              recv_2d_DP(:,:,:),    &
                              recvcount,            &
                              datatype,             &
                              PRC_masterrank,       &
                              PRC_LOCAL_COMM_WORLD, &
                              ierr                  )

             p = 1
             do py = 1, nprocs_y_out
             do px = 1, nprocs_x_out
                ! rearrangement of data array
                do j = 1, D2
                do i = 1, D1
                   dinfo_all%VAR_2d( i + (px-1) * D1, &
                                     j + (py-1) * D2  ) = recv_2d_DP(i,j,p)
                enddo
                enddo
                p = p + 1
             enddo
             enddo

             deallocate( send_2d_DP )
             deallocate( recv_2d_DP )

          endselect

       case( 3 )

          D1 = size( dinfo%VAR_3d(:,:,:), 1 )
          D2 = size( dinfo%VAR_3d(:,:,:), 2 )
          D3 = size( dinfo%VAR_3d(:,:,:), 3 )

          sendcount = D1 * D2 * D3
          recvcount = D1 * D2 * D3

          if ( dinfo_all%transpose ) then

             allocate( dinfo_all%VAR_3d( D1, D2*nprocs_x_out, D3*nprocs_y_out ) )

             dinfo_all%dim_size(1) = D1
             dinfo_all%dim_size(2) = D2 * nprocs_x_out
             dinfo_all%dim_size(3) = D3 * nprocs_y_out

             select case( dinfo_all%datatype )
             case( FILE_REAL4 )

                datatype = MPI_REAL

                allocate( send_3d_SP(D1,D2,D3       ) )
                allocate( recv_3d_SP(D1,D2,D3,nprocs) )

                send_3d_SP(:,:,:) = real( dinfo%VAR_3d(:,:,:), kind=SP )

                call MPI_GATHER( send_3d_SP(:,:,:),    &
                                 sendcount,            &
                                 datatype,             &
                                 recv_3d_SP(:,:,:,:),  &
                                 recvcount,            &
                                 datatype,             &
                                 PRC_masterrank,       &
                                 PRC_LOCAL_COMM_WORLD, &
                                 ierr                  )

                p = 1
                do py = 1, nprocs_y_out
                do px = 1, nprocs_x_out
                   ! rearrangement of data array
                   do j = 1, D3
                   do i = 1, D2
                   do k = 1, D1
                      dinfo_all%VAR_3d( k,               &
                                        i + (px-1) * D2, &
                                        j + (py-1) * D3  ) = recv_3d_SP(k,i,j,p)
                   enddo
                   enddo
                   enddo
                   p = p + 1
                enddo
                enddo

                deallocate( send_3d_SP )
                deallocate( recv_3d_SP )

             case( FILE_REAL8 )

                datatype = MPI_DOUBLE_PRECISION

                allocate( send_3d_DP(D1,D2,D3       ) )
                allocate( recv_3d_DP(D1,D2,D3,nprocs) )

                send_3d_DP(:,:,:) = real( dinfo%VAR_3d(:,:,:), kind=DP )

                call MPI_GATHER( send_3d_DP(:,:,:),    &
                                 sendcount,            &
                                 datatype,             &
                                 recv_3d_DP(:,:,:,:),  &
                                 recvcount,            &
                                 datatype,             &
                                 PRC_masterrank,       &
                                 PRC_LOCAL_COMM_WORLD, &
                                 ierr                  )

                p = 1
                do py = 1, nprocs_y_out
                do px = 1, nprocs_x_out
                   ! rearrangement of data array
                   do j = 1, D3
                   do i = 1, D2
                   do k = 1, D1
                      dinfo_all%VAR_3d( k,               &
                                        i + (px-1) * D2, &
                                        j + (py-1) * D3  ) = recv_3d_DP(k,i,j,p)
                   enddo
                   enddo
                   enddo
                   p = p + 1
                enddo
                enddo

                deallocate( send_3d_DP )
                deallocate( recv_3d_DP )

             endselect

          else

             allocate( dinfo_all%VAR_3d( D1*nprocs_x_out, D2*nprocs_y_out, D3 ) )

             dinfo_all%dim_size(1) = D1 * nprocs_x_out
             dinfo_all%dim_size(2) = D2 * nprocs_y_out
             dinfo_all%dim_size(3) = D3

             select case( dinfo_all%datatype )
             case( FILE_REAL4 )

                datatype = MPI_REAL

                allocate( send_3d_SP(D1,D2,D3       ) )
                allocate( recv_3d_SP(D1,D2,D3,nprocs) )

                send_3d_SP(:,:,:) = real( dinfo%VAR_3d(:,:,:), kind=SP )

                call MPI_GATHER( send_3d_SP(:,:,:),    &
                                 sendcount,            &
                                 datatype,             &
                                 recv_3d_SP(:,:,:,:),  &
                                 recvcount,            &
                                 datatype,             &
                                 PRC_masterrank,       &
                                 PRC_LOCAL_COMM_WORLD, &
                                 ierr                  )

                p = 1
                do py = 1, nprocs_y_out
                do px = 1, nprocs_x_out
                   ! rearrangement of data array
                   do k = 1, D3
                   do j = 1, D2
                   do i = 1, D1
                      dinfo_all%VAR_3d( i + (px-1) * D1, &
                                        j + (py-1) * D2, &
                                        k                ) = recv_3d_SP(i,j,k,p)
                   enddo
                   enddo
                   enddo
                   p = p + 1
                enddo
                enddo

                deallocate( send_3d_SP )
                deallocate( recv_3d_SP )

             case( FILE_REAL8 )

                datatype = MPI_DOUBLE_PRECISION

                allocate( send_3d_DP(D1,D2,D3       ) )
                allocate( recv_3d_DP(D1,D2,D3,nprocs) )

                send_3d_DP(:,:,:) = real( dinfo%VAR_3d(:,:,:), kind=DP )

                call MPI_GATHER( send_3d_DP(:,:,:),    &
                                 sendcount,            &
                                 datatype,             &
                                 recv_3d_DP(:,:,:,:),  &
                                 recvcount,            &
                                 datatype,             &
                                 PRC_masterrank,       &
                                 PRC_LOCAL_COMM_WORLD, &
                                 ierr                  )

                p = 1
                do py = 1, nprocs_y_out
                do px = 1, nprocs_x_out
                   ! rearrangement of data array
                   do k = 1, D3
                   do j = 1, D2
                   do i = 1, D1
                      dinfo_all%VAR_3d( i + (px-1) * D1, &
                                        j + (py-1) * D2, &
                                        k                ) = recv_3d_DP(i,j,k,p)
                   enddo
                   enddo
                   enddo
                   p = p + 1
                enddo
                enddo

                deallocate( send_3d_DP )
                deallocate( recv_3d_DP )

             endselect

          endif

       endselect

    else
       ! copy variable information
       dinfo_all = dinfo

    endif

    return
  end subroutine SNO_comm_globalvars

end module mod_sno_comm
