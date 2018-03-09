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
  use scale_comm_cartesC_nest, only: &
       PARENT_KMAX,     &
       PARENT_IMAX,     &
       PARENT_JMAX,     &
       PARENT_LKMAX,    &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y => COMM_CARTESC_NEST_TILE_NUM_Y, &
       NEST_TILE_ID    => COMM_CARTESC_NEST_TILE_ID

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtmosSetupSCALE
  public :: ParentAtmosOpenSCALE
  public :: ParentAtmosInputSCALE
  public :: ParentLandSetupSCALE
  public :: ParentLandInputSCALE
  public :: ParentOceanSetupSCALE
  public :: ParentOceanOpenSCALE
  public :: ParentOceanInputSCALE

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
  !> Atmos Setup
  subroutine ParentAtmosSetupSCALE( &
       dims )
    implicit none

    integer,          intent(out) :: dims(6)

    integer :: i
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atmos Input File Type: SCALE-RM'
    ! full level
    dims(1) = PARENT_KMAX(handle)
    dims(2) = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    dims(3) = PARENT_JMAX(handle) * NEST_TILE_NUM_Y
    ! half level
    dims(4) = dims(1)
    dims(5) = dims(2)
    dims(6) = dims(3)

    allocate( read2D ( PARENT_IMAX(handle), PARENT_JMAX(handle) ) )
    allocate( read3D ( PARENT_IMAX(handle), PARENT_JMAX(handle), dims(1) ) )

    return
  end subroutine ParentAtmosSetupSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosOpenSCALE( &
       lon_org,      &
       lat_org,      &
       cz_org,       &
       basename_org, &
       dims          )
    use scale_const, only: &
         D2R => CONST_D2R
    use scale_file, only: &
         FILE_open, &
         FILE_read
    implicit none
    real(RP), intent(out) :: lon_org(:,:)
    real(RP), intent(out) :: lat_org(:,:)
    real(RP), intent(out) :: cz_org (:,:,:)
    character(len=*), intent(in)  :: basename_org
    integer,  intent(in)  :: dims(6)

    integer :: rank
    integer :: xloc, yloc
    integer :: xs, xe, ys, ye

    integer :: fid
    integer :: i, k

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[AtmosOpenSCALE]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX(handle) * (xloc-1) + 1
       xe = PARENT_IMAX(handle) * xloc
       ys = PARENT_JMAX(handle) * (yloc-1) + 1
       ye = PARENT_JMAX(handle) * yloc

       call FILE_open( BASENAME_ORG,                  & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       call FILE_read( fid, "lon", read2D(:,:) )
       lon_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FILE_read( fid, "lat", read2D(:,:) )
       lat_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FILE_read( fid, "height", read3D(:,:,:) )
       do k = 1, dims(1)
          cz_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FILE_read( fid, "topo", read2D(:,:) )
       cz_org(2,xs:xe,ys:ye)  = read2D(:,:)

    end do

    cz_org(1,:,:)  = 0.0_RP
    return
  end subroutine ParentAtmosOpenSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosInputSCALE( &
       velz_org,      &
       velx_org,      &
       vely_org,      &
       pres_org,      &
       dens_org,      &
       pott_org,      &
       qv_org,        &
       qhyd_org,      &
       qnum_org,      &
       qtrc_org,      &
       cz_org,        &
       basename_org,  &
       same_mptype,   &
       dims,          &
       it             ) ! (in)
    use scale_const, only: &
       P00 => CONST_PRE00, &
       CPdry => CONST_CPdry, &
       Rdry => CONST_Rdry, &
       GRAV => CONST_GRAV, &
       LAPS => CONST_LAPS
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_file, only: &
       FILE_open, &
       FILE_read
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_atmos_thermodyn, only: &
       THERMODYN_specific_heat  => ATMOS_THERMODYN_specific_heat, &
       THERMODYN_rhot2temp_pres => ATMOS_THERMODYN_rhot2temp_pres
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       HYD_NAME, &
       NUM_NAME
    use scale_atmos_grid_cartesC_metric, only: &
       rotc => ATMOS_GRID_CARTESC_METRIC_ROTC
    use scale_topography, only: &
       topo => TOPO_Zsfc
    implicit none

    real(RP),         intent(out) :: velz_org(:,:,:)
    real(RP),         intent(out) :: velx_org(:,:,:)
    real(RP),         intent(out) :: vely_org(:,:,:)
    real(RP),         intent(out) :: pres_org(:,:,:)
    real(RP),         intent(out) :: dens_org(:,:,:)
    real(RP),         intent(out) :: pott_org(:,:,:)
    real(RP),         intent(out) :: qv_org  (:,:,:)
    real(RP),         intent(out) :: qhyd_org(:,:,:,:)
    real(RP),         intent(out) :: qnum_org(:,:,:,:)
    real(RP),         intent(out) :: qtrc_org(:,:,:,:)
    real(RP),         intent(in)  :: cz_org(:,:,:)
    character(len=*), intent(in)  :: basename_org
    logical,          intent(in)  :: same_mptype
    integer,          intent(in)  :: dims(6)
    integer,          intent(in)  :: it

    ! work

    real(RP) :: momz_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: momx_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: momy_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: rhot_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: tsfc_org(          dims(2),dims(3))
    real(RP) :: Qdry, Rtot, CVtot, CPtot
    real(RP) :: temp_org
    real(RP) :: dz

    integer :: xs, xe
    integer :: ys, ye
    integer :: xloc, yloc
    integer :: rank

    integer :: fid
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX(handle) * (xloc-1) + 1
       xe = PARENT_IMAX(handle) * xloc
       ys = PARENT_JMAX(handle) * (yloc-1) + 1
       ye = PARENT_JMAX(handle) * yloc

       call FILE_open( BASENAME_ORG,                  & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       call FILE_read( fid, "T2", read2D(:,:), step=it )
!OCL XFILL
       tsfc_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_read( fid, "MSLP", read2D(:,:), step=it )
!OCL XFILL
       pres_org(1,xs:xe,ys:ye) = read2D(:,:)

       call FILE_read( fid, "DENS", read3D(:,:,:), step=it )
!OCL XFILL
       do k = 1, dims(1)
          dens_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FILE_read( fid, "MOMZ", read3D(:,:,:), step=it )
!OCL XFILL
       do k = 1, dims(1)
          momz_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FILE_read( fid, "MOMX", read3D(:,:,:), step=it )
!OCL XFILL
       do k = 1, dims(1)
          momx_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FILE_read( fid, "MOMY", read3D(:,:,:), step=it )
!OCL XFILL
       do k = 1, dims(1)
          momy_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FILE_read( fid, "RHOT", read3D(:,:,:), step=it )
!OCL XFILL
       do k = 1, dims(1)
          rhot_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
       end do

!OCL XFILL
       do iq = 1, N_HYD
          qhyd_org(:,xs:xe,ys:ye,iq) = 0.0_RP
          qnum_org(:,xs:xe,ys:ye,iq) = 0.0_RP
       end do

       if( same_mptype ) then

          do iq = QS_MP, QE_MP
             call FILE_read( fid, TRACER_NAME(iq), read3D(:,:,:), step=it )
             do k = 1, dims(1)
                qtrc_org(k+2,xs:xe,ys:ye,iq) = read3D(:,:,k)
             end do
!OCL XFILL
             qtrc_org(2,xs:xe,ys:ye,iq) = qtrc_org(3,xs:xe,ys:ye,iq)
!OCL XFILL
             qtrc_org(1,xs:xe,ys:ye,iq) = qtrc_org(3,xs:xe,ys:ye,iq)
          enddo

       else 

          call FILE_read( fid, "QV", read3D(:,:,:), step=it, allow_missing=.true. )
!OCL XFILL
          do k = 1, dims(1)
             qv_org(k+2,xs:xe,ys:ye) = read3D(:,:,k)
          end do
!OCL XFILL
          qv_org(2,xs:xe,ys:ye) = qv_org(3,xs:xe,ys:ye)
!OCL XFILL
          qv_org(1,xs:xe,ys:ye) = qv_org(3,xs:xe,ys:ye)

          do iq = 1, N_HYD
             call FILE_read( fid, HYD_NAME(iq), read3D(:,:,:), step=it, allow_missing=.true. )
!OCL XFILL
             do k = 1, dims(1)
                qhyd_org(k+2,xs:xe,ys:ye,iq) = read3D(:,:,k)
             end do
!OCL XFILL
             qhyd_org(2,xs:xe,ys:ye,iq) = qhyd_org(3,xs:xe,ys:ye,iq)
!OCL XFILL
             qhyd_org(1,xs:xe,ys:ye,iq) = qhyd_org(3,xs:xe,ys:ye,iq)

             ! number density
             call FILE_read( fid, NUM_NAME(iq), read3D(:,:,:), step=it, allow_missing=.true. )
!OCL XFILL
             do k = 1, dims(1)
                qnum_org(k+2,xs:xe,ys:ye,iq) = read3D(:,:,k)
             end do
!OCL XFILL
             qnum_org(2,xs:xe,ys:ye,iq) = qnum_org(3,xs:xe,ys:ye,iq)
!OCL XFILL
             qnum_org(1,xs:xe,ys:ye,iq) = qnum_org(3,xs:xe,ys:ye,iq)
          end do

       endif

!       call FILE_read( fid, "Q2", read2D(:,:), step=it )
!       qv_org(2,xs:xe,ys:ye) = read2D(:,:)

    end do

    ! convert from momentum to velocity
!OCL XFILL
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 4, dims(1)+2
       velz_org(k,i,j) = ( momz_org(k-1,i,j) + momz_org(k,i,j) ) / dens_org(k,i,j) * 0.5_RP
    end do
    end do
    end do
!OCL XFILL
    do j = 1, dims(3)
    do i = 1, dims(2)
       velz_org(1:3    ,i,j) = 0.0_RP
       velz_org(dims(1)+2,i,j) = 0.0_RP
    end do
    end do

    ! convert from momentum to velocity
!OCL XFILL
    do j = 1, dims(3)
    do i = 2, dims(2)
    do k = 3, dims(1)+2
       velx_org(k,i,j) = ( momx_org(k,i-1,j) + momx_org(k,i,j) ) / dens_org(k,i,j) * 0.5_RP
    end do
    end do
    end do
!OCL XFILL
    do j = 1, dims(3)
    do k = 3, dims(1)+2
       velx_org(k,1,j) = momx_org(k,1,j) / dens_org(k,1,j)
    end do
    end do
!OCL XFILL
    velx_org(1:2,:,:) = 0.0_RP

    ! convert from momentum to velocity
!OCL XFILL
    do j = 2, dims(3)
    do i = 1, dims(2)
    do k = 3, dims(1)+2
       vely_org(k,i,j) = ( momy_org(k,i,j-1) + momy_org(k,i,j) ) / dens_org(k,i,j) * 0.5_RP
    end do
    end do
    end do
!OCL XFILL
    do i = 1, dims(2)
    do k = 3, dims(1)+2
       vely_org(k,i,1) = momy_org(k,i,1) / dens_org(k,i,1)
    end do
    end do
!OCL XFILL
    vely_org(1:2,:,:) = 0.0_RP


    !!! must be rotate !!!


    ! diagnose temp and pres
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 3, dims(1)+2
       call THERMODYN_specific_heat( QA, &
                                     qtrc_org(k,i,j,:), &
                                     TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! [IN]
                                     Qdry, Rtot, CVtot, CPtot                                 ) ! [OUT]
       call THERMODYN_rhot2temp_pres( dens_org(k,i,j), rhot_org(k,i,j), Rtot, CVtot, CPtot, &
                                      temp_org, pres_org(k,i,j)                             )
    end do
    end do
    end do
!OCL XFILL
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 3, dims(1)+2
       pott_org(k,i,j) = rhot_org(k,i,j) / dens_org(k,i,j)
    end do
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)
       dz = cz_org(3,i,j) - cz_org(2,i,j)

       dens_org(2,i,j) = ( pres_org(3,i,j) + GRAV * dens_org(3,i,j) * dz * 0.5_RP ) &
                       / ( Rdry * tsfc_org(i,j) - GRAV * dz * 0.5_RP )
       pres_org(2,i,j) = dens_org(2,i,j) * Rdry * tsfc_org(i,j)
       pott_org(2,i,j) = tsfc_org(i,j) * ( P00 / pres_org(2,i,j) )**(Rdry/CPdry)

       temp_org = tsfc_org(i,j) + LAPS * cz_org(2,i,j)
       pott_org(1,i,j) = temp_org * ( P00 / pres_org(1,i,j) )**(Rdry/CPdry)
       dens_org(1,i,j) = pres_org(1,i,j) / ( Rdry * temp_org )
    end do
    end do

    return
  end subroutine ParentAtmosInputSCALE

  !-----------------------------------------------------------------------------
  !> Land Setup
  subroutine ParentLandSetupSCALE( &
       ldims )
    implicit none

    integer, intent(out) :: ldims(3)

    integer :: i
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Land Input File Type: SCALE-RM'
    ldims(1) = PARENT_LKMAX(handle)
    ldims(2) = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    ldims(3) = PARENT_JMAX(handle) * NEST_TILE_NUM_Y

    if ( .not. allocated(read2D) ) then
       allocate( read2D ( PARENT_IMAX(handle), PARENT_JMAX(handle) ) )
    end if
    allocate( read3DL( PARENT_IMAX(handle), PARENT_JMAX(handle), ldims(1) ) )

    return
  end subroutine ParentLandSetupSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentLandInputSCALE( &
      tg_org,             &
      strg_org,           &
      lst_org,            &
      ust_org,            &
      albg_org,           &
      topo_org,           &
      lmask_org,          &
      llon_org,           &
      llat_org,           &
      lz_org,             &
      basename_land,      &
      ldims,              &
      use_file_landwater, &
      it                  )
    use scale_file, only: &
         FILE_open, &
         FILE_read
    use scale_const, only: &
         D2R => CONST_D2R
    implicit none

    real(RP), intent(out) :: tg_org(:,:,:)
    real(RP), intent(out) :: strg_org(:,:,:)
    real(RP), intent(out) :: lst_org(:,:)
    real(RP), intent(out) :: ust_org(:,:)
    real(RP), intent(out) :: albg_org(:,:,:)
    real(RP), intent(out) :: topo_org(:,:)
    real(RP), intent(out) :: lmask_org(:,:)
    real(RP), intent(out) :: llon_org(:,:)
    real(RP), intent(out) :: llat_org(:,:)
    real(RP), intent(out) :: lz_org(:)

    character(len=*), intent(in) :: basename_land
    integer,          intent(in) :: ldims(3)
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    integer,          intent(in) :: it

    integer :: rank

    integer :: fid
    integer :: k, i, j, n
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[LandInputSCALE]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX(handle) * (xloc-1) + 1
       xe = PARENT_IMAX(handle) * xloc
       ys = PARENT_JMAX(handle) * (yloc-1) + 1
       ye = PARENT_JMAX(handle) * yloc

       call FILE_open( BASENAME_land,                 & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       call FILE_read( fid, "LAND_TEMP",  read3DL(:,:,:), step=it )
       do k = 1, ldims(1)
         tg_org(k,xs:xe,ys:ye) = read3DL(:,:,k)
       end do

       if( use_file_landwater )then
          call FILE_read( fid, "LAND_WATER", read3DL(:,:,:), step=it )
          do k = 1, ldims(1)
             strg_org(k,xs:xe,ys:ye) = read3DL(:,:,k)
          end do
       endif

       call FILE_read( fid, "lon", read2D(:,:) )
       llon_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FILE_read( fid, "lat", read2D(:,:) )
       llat_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FILE_read( fid, "LAND_SFC_TEMP", read2D(:,:), step=it )
       lst_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_read( fid, "URBAN_SFC_TEMP", read2D(:,:), step=it )
       ust_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_read( fid, "LAND_ALB_LW", read2D(:,:), step=it )
       albg_org(xs:xe,ys:ye,1) = read2D(:,:)

       call FILE_read( fid, "LAND_ALB_SW", read2D(:,:), step=it )
       albg_org(xs:xe,ys:ye,2) = read2D(:,:)

       call FILE_read( fid, "topo", read2D(:,:) )
       topo_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_read( fid, "lsmask", read2D(:,:) )
       lmask_org(xs:xe,ys:ye) = read2D(:,:)

    end do

    call FILE_read( fid, "lz", lz_org(:) )

    return
  end subroutine ParentLandInputSCALE

  !-----------------------------------------------------------------------------
  !> Ocean Setup
  subroutine ParentOceanSetupSCALE( &
       odims )
    implicit none

    integer,          intent(out) :: odims(2)

    integer :: i
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Ocean Input File Type: SCALE-RM'
    odims(1) = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    odims(2) = PARENT_JMAX(handle) * NEST_TILE_NUM_Y

    if ( .not. allocated(read2D) ) then
       allocate( read2D ( PARENT_IMAX(handle), PARENT_JMAX(handle) ) )
    end if

    return
  end subroutine ParentOceanSetupSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentOceanOpenSCALE( &
       olon_org,      &
       olat_org,      &
       omask_org, &
       basename_ocean, &
       odims          )
    use scale_const, only: &
         D2R => CONST_D2R
    use scale_file, only: &
         FILE_open, &
         FILE_read
    implicit none
    real(RP), intent(out) :: olon_org (:,:)
    real(RP), intent(out) :: olat_org (:,:)
    real(RP), intent(out) :: omask_org(:,:)
    character(len=*), intent(in)  :: basename_ocean
    integer,  intent(in)  :: odims(2)

    integer :: rank
    integer :: xloc, yloc
    integer :: xs, xe, ys, ye

    integer :: fid
    integer :: i, k

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OceanOpenSCALE]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX(handle) * (xloc-1) + 1
       xe = PARENT_IMAX(handle) * xloc
       ys = PARENT_JMAX(handle) * (yloc-1) + 1
       ye = PARENT_JMAX(handle) * yloc

       call FILE_open( BASENAME_ocean,                & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       call FILE_read( fid, "lon", read2D(:,:) )
       olon_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FILE_read( fid, "lat", read2D(:,:) )
       olat_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FILE_read( fid, "lsmask", read2D(:,:) )
       omask_org(xs:xe,ys:ye)  = read2D(:,:)

    end do

    return
  end subroutine ParentOceanOpenSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentOceanInputSCALE( &
      tw_org,         &
      sst_org,        &
      albw_org,       &
      z0w_org,        &
      omask_org,      &
      basename_ocean, &
      odims,          &
      it              )
    use scale_file, only: &
         FILE_open, &
         FILE_read
    implicit none

    real(RP), intent(out) :: tw_org(:,:)
    real(RP), intent(out) :: sst_org(:,:)
    real(RP), intent(out) :: albw_org(:,:,:)
    real(RP), intent(out) :: z0w_org(:,:)
    real(RP), intent(out) :: omask_org(:,:)

    character(len=*), intent(in) :: basename_ocean
    integer,          intent(in) :: odims(2)
    integer,          intent(in) :: it

    integer :: rank

    integer :: fid
    integer :: i, j, n
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OceanInputSCALE]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX(handle) * (xloc-1) + 1
       xe = PARENT_IMAX(handle) * xloc
       ys = PARENT_JMAX(handle) * (yloc-1) + 1
       ye = PARENT_JMAX(handle) * yloc

       call FILE_open( BASENAME_ocean,                & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       call FILE_read( fid, "OCEAN_TEMP", read2D(:,:), step=it )
       tw_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_read( fid, "OCEAN_SFC_TEMP", read2D(:,:), step=it )
       sst_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_read( fid, "OCEAN_ALB_LW", read2D(:,:), step=it )
       albw_org(xs:xe,ys:ye,1) = read2D(:,:)

       call FILE_read( fid, "OCEAN_ALB_SW", read2D(:,:), step=it )
       albw_org(xs:xe,ys:ye,2) = read2D(:,:)

       call FILE_read( fid, "OCEAN_SFC_Z0M", read2D(:,:), step=it )
       z0w_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_read( fid, "lsmask", read2D(:,:) )
       omask_org(xs:xe,ys:ye) = read2D(:,:)

    end do

    return
  end subroutine ParentOceanInputSCALE

end module mod_realinput_scale
