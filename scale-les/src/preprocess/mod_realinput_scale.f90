!-------------------------------------------------------------------------------
!> module REAL input SCALE
!!
!! @par Description
!!          read data from SCALE file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-05-24 (S.Nishizawa)   [new] split from mod_realinput.f90
!!
!<
module mod_realinput_scale
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_tracer
  use scale_process, only: &
       myrank => PRC_myrank,  &
       PRC_MPIstop
  use scale_grid_nest, only: &
       PARENT_KMAX,     &
       PARENT_IMAX,     &
       PARENT_JMAX,     &
       PARENT_LKMAX,    &
       NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y, &
       NEST_TILE_ID

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtomSetupSCALE
  public :: ParentAtomOpenSCALE
  public :: ParentAtomInputSCALE
  public :: ParentSurfaceInputSCALE

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
  integer, parameter    :: handle = 1

  real(RP), allocatable :: read2D(:,:)
  real(RP), allocatable :: read3D(:,:,:)
  real(RP), allocatable :: read3DL(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ParentAtomSetupSCALE( &
       dims, &
       timelen )
    implicit none

    integer,          intent(out) :: dims(11)
    integer,          intent(out) :: timelen

    integer :: i
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: SCALE-LES'
    ! full level
    dims(1) = PARENT_KMAX(handle)
    dims(2) = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    dims(3) = PARENT_JMAX(handle) * NEST_TILE_NUM_Y
    ! half level
    dims(4) = dims(1)
    dims(5) = dims(2)
    dims(6) = dims(3)
    ! land
    dims(7) = PARENT_LKMAX(handle)
    dims(8) = dims(2)
    dims(9) = dims(3)
    ! sst
    dims(10) = dims(2)
    dims(11) = dims(3)

    timelen = 1 !temtative approach

    allocate( read2D ( PARENT_IMAX(handle), PARENT_JMAX(handle) ) )
    allocate( read3D ( PARENT_IMAX(handle), PARENT_JMAX(handle), dims(1) ) )
    allocate( read3DL( PARENT_IMAX(handle), PARENT_JMAX(handle), dims(7) ) )

    return
  end subroutine ParentAtomSetupSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentAtomOpenSCALE( &
       lon_org,      &
       lat_org,      &
       cz_org,       &
       basename_org, &
       dims          )
    use scale_const, only: &
         D2R => CONST_D2R
    use gtool_file, only: &
         FileRead
    implicit none
    real(RP), intent(out) :: lon_org(:,:)
    real(RP), intent(out) :: lat_org(:,:)
    real(RP), intent(out) :: cz_org (:,:,:)
    character(len=*), intent(in)  :: basename_org
    integer,  intent(in)  :: dims(11)

    integer :: rank
    integer :: xloc, yloc
    integer :: xs, xe, ys, ye

    integer :: i, k

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OpenSCALE]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX(handle) * (xloc-1) + 1
       xe = PARENT_IMAX(handle) * xloc
       ys = PARENT_JMAX(handle) * (yloc-1) + 1
       ye = PARENT_JMAX(handle) * yloc

       call FileRead( read2D(:,:),   BASENAME_ORG, "lon",        1, rank )
       lon_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FileRead( read2D(:,:),   BASENAME_ORG, "lat",        1, rank )
       lat_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FileRead( read3D(:,:,:), BASENAME_ORG, "height",     1, rank )
       do k = 1, dims(1)
          cz_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FileRead( read2D(:,:),   BASENAME_ORG, "topo",       1, rank )
       cz_org(2,xs:xe,ys:ye)  = read2D(:,:)

    end do

    cz_org(1,:,:)  = 0.0_RP
    return
  end subroutine ParentAtomOpenSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentAtomInputSCALE( &
       velz_org,     &
       velx_org,     &
       vely_org,     &
       pres_org,     &
       dens_org,     &
       pott_org,     &
       qtrc_org,     &
       basename_org, &
       dims,         &
       it            ) ! (in)
    use scale_const, only: &
       P00 => CONST_PRE00, &
       CPdry => CONST_CPdry, &
       Rdry => CONST_Rdry
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use gtool_file, only: &
       FileRead
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
       THERMODYN_pott      => ATMOS_THERMODYN_pott
    use scale_gridtrans, only: &
       rotc => GTRANS_ROTC
    use scale_topography, only: &
       topo => TOPO_Zsfc
    use scale_interpolation_nest, only: &
       INTRPNEST_interp_fact_llz,  &
       INTRPNEST_interp_2d,        &
       INTRPNEST_interp_3d
    implicit none

    real(RP),         intent(out) :: velz_org(:,:,:)
    real(RP),         intent(out) :: velx_org(:,:,:)
    real(RP),         intent(out) :: vely_org(:,:,:)
    real(RP),         intent(out) :: pres_org(:,:,:)
    real(RP),         intent(out) :: dens_org(:,:,:)
    real(RP),         intent(out) :: pott_org(:,:,:)
    real(RP),         intent(out) :: qtrc_org(:,:,:,:)
    character(len=*), intent(in)  :: basename_org
    integer,          intent(in)  :: dims(7)
    integer,          intent(in)  :: it

    ! work

    real(RP) :: momz_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: momx_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: momy_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: rhot_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: tsfc_org(          dims(2),dims(3))
    real(RP) :: temp_org

    integer :: xs, xe
    integer :: ys, ye
    integer :: xloc, yloc
    integer :: rank

    integer :: k, i, j, iq
    logical :: lack_of_val
    !---------------------------------------------------------------------------


    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputSCALE]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX(handle) * (xloc-1) + 1
       xe = PARENT_IMAX(handle) * xloc
       ys = PARENT_JMAX(handle) * (yloc-1) + 1
       ye = PARENT_JMAX(handle) * yloc

       call FileRead( read2D(:,:), BASENAME_ORG, "T2", it, rank )
       tsfc_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "MSLP", it, rank )
       pres_org(1,xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read3D(:,:,:), BASENAME_ORG, "DENS", it, rank )
       do k = 1, dims(1)
          dens_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FileRead( read3D(:,:,:), BASENAME_ORG, "MOMZ", it, rank )
       do k = 1, dims(1)
          momz_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FileRead( read3D(:,:,:), BASENAME_ORG, "MOMX", it, rank )
       do k = 1, dims(1)
          momx_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FileRead( read3D(:,:,:), BASENAME_ORG, "MOMY", it, rank )
       do k = 1, dims(1)
          momy_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FileRead( read3D(:,:,:), BASENAME_ORG, "RHOT", it, rank )
       do k = 1, dims(1)
          rhot_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       do iq = 1, QA
          call FileRead( read3D(:,:,:), BASENAME_ORG, AQ_NAME(iq), it, rank )
          do k = 1, dims(1)
             qtrc_org(k+2,xs:xe,ys:ye,iq) = read3D(:,:,k)
          end do
          qtrc_org(2,xs:xe,ys:ye,iq) = qtrc_org(3,xs:xe,ys:ye,iq)
       end do

!       call FileRead( read2D(:,:), BASENAME_ORG, "Q2", it, rank )
!       qtrc_org(2,xs:xe,ys:ye,I_QV) = read2D(:,:)

    end do

    ! convert from momentum to velocity
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 4, dims(1)+2
       velz_org(k,i,j) = ( momz_org(k-1,i,j) + momz_org(k,i,j) ) / dens_org(k,i,j) * 0.5_RP
    end do
    end do
    end do
    do j = 1, dims(3)
    do i = 1, dims(2)
       velz_org(1:3    ,i,j) = 0.0_RP
       velz_org(dims(1)+2,i,j) = 0.0_RP
    end do
    end do

    ! convert from momentum to velocity
    do j = 1, dims(3)
    do i = 2, dims(2)
    do k = 3, dims(1)+2
       velx_org(k,i,j) = ( momx_org(k,i-1,j) + momx_org(k,i,j) ) / dens_org(k,i,j) * 0.5_RP
    end do
    end do
    end do
    do j = 1, dims(3)
    do k = 3, dims(1)+2
       velx_org(k,1,j) = momx_org(k,1,j) / dens_org(k,1,j)
    end do
    end do
    velx_org(1:2,:,:) = 0.0_RP

    ! convert from momentum to velocity
    do j = 2, dims(3)
    do i = 1, dims(2)
    do k = 3, dims(1)+2
       vely_org(k,i,j) = ( momy_org(k,i,j-1) + momy_org(k,i,j) ) / dens_org(k,i,j) * 0.5_RP
    end do
    end do
    end do
    do i = 1, dims(2)
    do k = 3, dims(1)+2
       vely_org(k,i,1) = momy_org(k,i,1) / dens_org(k,i,1)
    end do
    end do
    vely_org(1:2,:,:) = 0.0_RP


    !!! must be rotate !!!


    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 3, dims(1)+2
          ! diagnose temp and pres
          call THERMODYN_temp_pres( temp_org,          & ! [OUT]
                                    pres_org(k,i,j),   & ! [OUT]
                                    dens_org(k,i,j),   & ! [IN]
                                    rhot_org(k,i,j),   & ! [IN]
                                    qtrc_org(k,i,j,:)  ) ! [IN]
          pott_org(k,i,j) = rhot_org(k,i,j) / dens_org(k,i,j)
       end do
       pott_org(1:2,i,j) = pott_org(3,i,j)
       pres_org(2,i,j) = P00 * ( tsfc_org(i,j) / pott_org(2,i,j) )**(CPdry/Rdry)
       dens_org(1,i,j) = P00 * ( P00/pres_org(1,i,j) )**(Rdry/CPdry-1.0_RP)
       dens_org(2,i,j) = pres_org(2,i,j) / ( tsfc_org(i,j) * Rdry )
       qtrc_org(1,i,j,:) = qtrc_org(2,i,j,:)
    end do
    end do

    return
  end subroutine ParentAtomInputSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentSurfaceInputSCALE( &
      tg_org,             &
      strg_org,           &
      tw_org,             &
      lst_org,            &
      ust_org,            &
      sst_org,            &
      albw_org,           &
      albg_org,           &
      z0w_org,            &
      lsmask_org,         &
      lz_org,             &
      basename_org,       &
      dims,               &
      use_file_landwater, &
      it                  )
    use gtool_file, only: &
         FileRead
    implicit none

    real(RP), intent(out) :: tg_org(:,:,:)
    real(RP), intent(out) :: strg_org(:,:,:)
    real(RP), intent(out) :: tw_org(:,:)
    real(RP), intent(out) :: lst_org(:,:)
    real(RP), intent(out) :: ust_org(:,:)
    real(RP), intent(out) :: sst_org(:,:)
    real(RP), intent(out) :: albw_org(:,:,:)
    real(RP), intent(out) :: albg_org(:,:,:)
    real(RP), intent(out) :: z0w_org(:,:)
    real(RP), intent(out) :: lsmask_org(:,:)
    real(RP), intent(out) :: lz_org(:)

    character(len=*), intent(in) :: basename_org
    integer,          intent(in) :: dims(11)
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    integer,          intent(in) :: it

    integer :: rank

    integer :: k, i, j, n
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputSCALE-Surface]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX(handle) * (xloc-1) + 1
       xe = PARENT_IMAX(handle) * xloc
       ys = PARENT_JMAX(handle) * (yloc-1) + 1
       ye = PARENT_JMAX(handle) * yloc

       call FileRead( read3DL(:,:,:), BASENAME_ORG, "LAND_TEMP",  it, rank )
       do k = 1, dims(7)
         tg_org(k,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       if( use_file_landwater )then
          call FileRead( read3DL(:,:,:), BASENAME_ORG, "LAND_WATER", it, rank )
          do k = 1, dims(7)
             strg_org(k,xs:xe,ys:ye) = read3D(:,:,k)
          end do
       endif

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_TEMP",     it, rank )
       tw_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "LAND_SFC_TEMP",  it, rank )
       lst_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "URBAN_SFC_TEMP", it, rank )
       ust_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_SFC_TEMP", it, rank )
       sst_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_ALB_LW",   it, rank )
       albw_org(xs:xe,ys:ye,1) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_ALB_SW",   it, rank )
       albw_org(xs:xe,ys:ye,2) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "LAND_ALB_LW",    it, rank )
       albg_org(xs:xe,ys:ye,1) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "LAND_ALB_SW",    it, rank )
       albg_org(xs:xe,ys:ye,2) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_SFC_Z0M",  it, rank )
       z0w_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "lsmask",         it, rank )
       lsmask_org(xs:xe,ys:ye) = read2D(:,:)

    end do

    call FileRead( lz_org(:),   BASENAME_ORG, "lz",  1, rank )

    return
  end subroutine ParentSurfaceInputSCALE

end module mod_realinput_scale
!-------------------------------------------------------------------------------
