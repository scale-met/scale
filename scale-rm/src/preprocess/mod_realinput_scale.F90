!-------------------------------------------------------------------------------
!> module REAL input SCALE
!!
!! @par Description
!!          read data from SCALE file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module mod_realinput_scale
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_tracer
  use scale_cpl_sfc_index
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
  integer,  private, parameter   :: handle = 1

  real(RP), private, allocatable :: read2D (:,:)
  real(RP), private, allocatable :: read3D (:,:,:)
  real(RP), private, allocatable :: read3DL(:,:,:)
  real(RP), private, allocatable :: read3DO(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtmosSetupSCALE( &
       dims )
    use scale_comm_cartesC_nest, only: &
       PARENT_KMAX,                                     &
       PARENT_IMAX,                                     &
       PARENT_JMAX,                                     &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y => COMM_CARTESC_NEST_TILE_NUM_Y
    implicit none

    integer, intent(out) :: dims(6)
    !---------------------------------------------------------------------------

    LOG_INFO("ParentAtmosSetupSCALE",*) 'Real Case/Atmos Input File Type: SCALE-RM'

    ! full level
    dims(1) = PARENT_KMAX(handle)
    dims(2) = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    dims(3) = PARENT_JMAX(handle) * NEST_TILE_NUM_Y
    ! half level
    dims(4) = dims(1)
    dims(5) = dims(2)
    dims(6) = dims(3)

    allocate( read2D(        PARENT_IMAX(handle),PARENT_JMAX(handle)) )
    allocate( read3D(dims(1),PARENT_IMAX(handle),PARENT_JMAX(handle)) )

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
    use scale_comm_cartesC_nest, only: &
       PARENT_IMAX,                                     &
       PARENT_JMAX,                                     &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_ID    => COMM_CARTESC_NEST_TILE_ID
    use scale_file, only: &
       FILE_open
    use scale_file_CARTESC, only: &
       FILE_CARTESC_read
    implicit none

    real(RP),         intent(out) :: lon_org(:,:)
    real(RP),         intent(out) :: lat_org(:,:)
    real(RP),         intent(out) :: cz_org (:,:,:)
    character(len=*), intent(in)  :: basename_org
    integer,          intent(in)  :: dims(6)

    integer :: rank
    integer :: xloc, yloc
    integer :: xs, xe, ys, ye

    logical :: existed

    integer :: fid
    integer :: i
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

       call FILE_CARTESC_read( fid, "lon", read2D(:,:) )
       lon_org(xs:xe,ys:ye) = read2D(:,:) * D2R

       call FILE_CARTESC_read( fid, "lat", read2D(:,:) )
       lat_org(xs:xe,ys:ye) = read2D(:,:) * D2R

       call FILE_CARTESC_read( fid, "height", read3D(:,:,:) )
       cz_org(3:dims(1)+2,xs:xe,ys:ye) = read3D(:,:,:)

       call FILE_CARTESC_read( fid, "topo", read2D(:,:), existed=existed )
       cz_org(2,xs:xe,ys:ye) = read2D(:,:)

    end do

    cz_org(1,:,:) = 0.0_RP

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
       qtrc_org,      &
       cz_org,        &
       basename_org,  &
       sfc_diagnoses, &
       same_mptype,   &
       dims,          &
       it             )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       P00   => CONST_PRE00, &
       CPdry => CONST_CPdry, &
       Rdry  => CONST_Rdry,  &
       GRAV  => CONST_GRAV,  &
       LAPS  => CONST_LAPS
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_comm_cartesC_nest, only: &
       PARENT_IMAX,                                     &
       PARENT_JMAX,                                     &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_ID    => COMM_CARTESC_NEST_TILE_ID
    use scale_file, only: &
       FILE_open
    use scale_file_CARTESC, only: &
       FILE_CARTESC_read
    use scale_atmos_grid_cartesC_metric, only: &
       rotc => ATMOS_GRID_CARTESC_METRIC_ROTC
    use scale_atmos_thermodyn, only: &
       THERMODYN_specific_heat  => ATMOS_THERMODYN_specific_heat, &
       THERMODYN_rhot2temp_pres => ATMOS_THERMODYN_rhot2temp_pres
    use scale_atmos_hydrometeor, only: &
       N_HYD,    &
       HYD_NAME, &
       NUM_NAME
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    implicit none

    real(RP),         intent(out) :: velz_org(:,:,:)
    real(RP),         intent(out) :: velx_org(:,:,:)
    real(RP),         intent(out) :: vely_org(:,:,:)
    real(RP),         intent(out) :: pres_org(:,:,:)
    real(RP),         intent(out) :: dens_org(:,:,:)
    real(RP),         intent(out) :: pott_org(:,:,:)
    real(RP),         intent(out) :: qv_org  (:,:,:)
    real(RP),         intent(out) :: qtrc_org(:,:,:,:)
    real(RP),         intent(in)  :: cz_org  (:,:,:)
    character(len=*), intent(in)  :: basename_org
    logical,          intent(in)  :: sfc_diagnoses
    logical,          intent(in)  :: same_mptype
    integer,          intent(in)  :: dims(6)
    integer,          intent(in)  :: it

    real(RP) :: momz_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: momx_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: momy_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: rhot_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: tsfc_org(          dims(2),dims(3))
    real(RP) :: qhyd_org(dims(1)+2,dims(2),dims(3),N_HYD)
    real(RP) :: qnum_org(dims(1)+2,dims(2),dims(3),N_HYD)
    real(RP) :: qdry    (dims(1)+2)
    real(RP) :: Rtot    (dims(1)+2)
    real(RP) :: CVtot   (dims(1)+2)
    real(RP) :: CPtot   (dims(1)+2)
    real(RP) :: temp_org(dims(1)+2)
    real(RP) :: dz

    integer  :: xs, xe
    integer  :: ys, ye
    integer  :: nx, ny
    integer  :: xloc, yloc
    integer  :: rank

    integer  :: fid
    logical  :: existed, existed_t2, existed_mslp
    integer  :: k, i, j, iq
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

       nx = xe - xs + 1
       ny = ye - ys + 1

       call FILE_open( BASENAME_ORG,                  & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       if ( sfc_diagnoses ) then
          call FILE_CARTESC_read( fid, "T2", read2D(:,:), step=it, existed=existed_t2 )
          if ( existed_t2 ) then
!OCL XFILL
             tsfc_org(xs:xe,ys:ye) = read2D(:,:)
          end if

          call FILE_CARTESC_read( fid, "MSLP", read2D(:,:), step=it, existed=existed_mslp )
          if ( existed_mslp ) then
!OCL XFILL
             pres_org(1,xs:xe,ys:ye) = read2D(:,:)
          end if
       end if

       call FILE_CARTESC_read( fid, "DENS", read3D(:,:,:), step=it )
!OCL XFILL
       dens_org(3:dims(1)+2,xs:xe,ys:ye) = read3D(:,:,:)

       call FILE_CARTESC_read( fid, "MOMZ", read3D(:,:,:), step=it )
!OCL XFILL
       momz_org(3:dims(1)+2,xs:xe,ys:ye) = read3D(:,:,:)

       call FILE_CARTESC_read( fid, "MOMX", read3D(:,:,:), step=it )
!OCL XFILL
       momx_org(3:dims(1)+2,xs:xe,ys:ye) = read3D(:,:,:)

       call FILE_CARTESC_read( fid, "MOMY", read3D(:,:,:), step=it )
!OCL XFILL
       momy_org(3:dims(1)+2,xs:xe,ys:ye) = read3D(:,:,:)


       call FILE_CARTESC_read( fid, "RHOT", read3D(:,:,:), step=it )
!OCL XFILL
       rhot_org(3:dims(1)+2,xs:xe,ys:ye) = read3D(:,:,:)


       if ( same_mptype ) then

          do iq = QS_MP, QE_MP
             call FILE_CARTESC_read( fid, TRACER_NAME(iq), read3D(:,:,:), step=it )
!OCL XFILL
             qtrc_org(3:dims(1)+2,xs:xe,ys:ye,iq) = read3D(:,:,:)
!OCL XFILL
             qtrc_org(2,xs:xe,ys:ye,iq) = qtrc_org(3,xs:xe,ys:ye,iq)
!OCL XFILL
             qtrc_org(1,xs:xe,ys:ye,iq) = qtrc_org(3,xs:xe,ys:ye,iq)
          enddo

       else

!OCL XFILL
          do iq = 1, N_HYD
             qhyd_org(:,xs:xe,ys:ye,iq) = 0.0_RP
             qnum_org(:,xs:xe,ys:ye,iq) = 0.0_RP
          end do

          call FILE_CARTESC_read( fid, "QV", read3D(:,:,:), step=it, existed=existed )
          if ( existed ) then
!OCL XFILL
             qv_org(3:dims(1)+2,xs:xe,ys:ye) = read3D(:,:,:)
!OCL XFILL
             qv_org(2,xs:xe,ys:ye) = qv_org(3,xs:xe,ys:ye)
!OCL XFILL
             qv_org(1,xs:xe,ys:ye) = qv_org(3,xs:xe,ys:ye)
          else
!OCL XFILL
             qv_org(:,:,:) = 0.0_RP
          end if

          do iq = 1, N_HYD
             call FILE_CARTESC_read( fid, HYD_NAME(iq), read3D(:,:,:), step=it, existed=existed )
             if ( existed ) then
!OCL XFILL
                qhyd_org(3:dims(1)+2,xs:xe,ys:ye,iq) = read3D(:,:,:)
!OCL XFILL
                qhyd_org(2,xs:xe,ys:ye,iq) = qhyd_org(3,xs:xe,ys:ye,iq)
!OCL XFILL
                qhyd_org(1,xs:xe,ys:ye,iq) = qhyd_org(3,xs:xe,ys:ye,iq)
             else
!OCL XFILL
                qhyd_org(:,:,:,iq) = 0.0_RP
             end if

             ! number density
             call FILE_CARTESC_read( fid, NUM_NAME(iq), read3D(:,:,:), step=it, existed=existed )
             if ( existed ) then
!OCL XFILL
                qnum_org(3:dims(1)+2,xs:xe,ys:ye,iq) = read3D(:,:,:)
!OCL XFILL
                qnum_org(2,xs:xe,ys:ye,iq) = qnum_org(3,xs:xe,ys:ye,iq)
!OCL XFILL
                qnum_org(1,xs:xe,ys:ye,iq) = qnum_org(3,xs:xe,ys:ye,iq)
             else
!OCL XFILL
                qnum_org(:,:,:,iq) = 0.0_RP
             end if
          end do

       endif

!       call FILE_CARTESC_read( fid, "Q2", read2D(:,:), step=it )
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
    if ( sfc_diagnoses ) then
!OCL XFILL
       velx_org(1:2,:,:) = 0.0_RP
    else
!OCL XFILL
       velx_org(1:2,:,:) = UNDEF
    end if

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
    if ( sfc_diagnoses ) then
!OCL XFILL
       vely_org(1:2,:,:) = 0.0_RP
    else
!OCL XFILL
       vely_org(1:2,:,:) = UNDEF
    end if


    !!! must be rotate !!!


    do iq = 1, size(qtrc_org,4)
       if ( iq >= QS_MP .and. iq <= QE_MP ) cycle
!OCL XFILL
       qtrc_org(:,:,:,iq) = 0.0_RP
    end do

    if ( QA_MP > 0 .AND. .NOT. same_mptype ) then
       call ATMOS_PHY_MP_driver_qhyd2qtrc( dims(1)+2, 1, dims(1)+2, dims(2), 1, dims(2), dims(3), 1, dims(3), &
                                           qv_org(:,:,:), qhyd_org(:,:,:,:), & ! [IN]
                                           qtrc_org(:,:,:,QS_MP:QE_MP),      & ! [OUT]
                                           QNUM=qnum_org(:,:,:,:)            ) ! [IN]
    end if

    !$omp parallel do default(none) &
    !$omp private(dz,temp_org,qdry,rtot,cvtot,cptot) &
    !$omp shared(dims,QA,UNDEF,LAPS,P00,Rdry,CPdry,GRAV,TRACER_MASS,TRACER_R,TRACER_CV,TRACER_CP, &
    !$omp        dens_org,rhot_org,qtrc_org,pres_org,pott_org,tsfc_org,cz_org, &
    !$omp        existed_t2,existed_mslp,sfc_diagnoses)
    do j = 1, dims(3)
    do i = 1, dims(2)

       ! diagnose temp and pres
       call THERMODYN_specific_heat( dims(1)+2, 3, dims(1)+2, QA,                             & ! [IN]
                                     qtrc_org(:,i,j,:),                                       & ! [IN]
                                     TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! [IN]
                                     qdry(:), Rtot(:), CVtot(:), CPtot(:)                     ) ! [OUT]

       call THERMODYN_rhot2temp_pres( dims(1), 1, dims(1),                                       & ! [IN]
                                      dens_org(3:dims(1)+2,i,j), rhot_org(3:dims(1)+2,i,j),      & ! [IN]
                                      Rtot(3:dims(1)+2), CVtot(3:dims(1)+2), CPtot(3:dims(1)+2), & ! [IN]
                                      temp_org(3:dims(1)+2), pres_org(3:dims(1)+2,i,j)           ) ! [OUT]

!OCL XFILL
       do k = 3, dims(1)+2
          pott_org(k,i,j) = rhot_org(k,i,j) / dens_org(k,i,j)
       end do

       if ( sfc_diagnoses ) then
          ! at the surface
          dz = cz_org(3,i,j) - cz_org(2,i,j)

          if ( .not. existed_t2 ) then
             tsfc_org(i,j) = temp_org(3) + LAPS * dz
          end if
          dens_org(2,i,j) = ( pres_org(3,i,j) + GRAV * dens_org(3,i,j) * dz * 0.5_RP ) &
                          / ( Rdry * tsfc_org(i,j) - GRAV * dz * 0.5_RP )
          pres_org(2,i,j) = dens_org(2,i,j) * Rdry * tsfc_org(i,j)
          pott_org(2,i,j) = tsfc_org(i,j) * ( P00 / pres_org(2,i,j) )**(Rdry/CPdry)

          ! at the sea-level
          temp_org(1) = tsfc_org(i,j) + LAPS * cz_org(2,i,j)
          if ( existed_mslp ) then
             pott_org(1,i,j) = temp_org(1) * ( P00 / pres_org(1,i,j) )**(Rdry/CPdry)
             dens_org(1,i,j) = pres_org(1,i,j) / ( Rdry * temp_org(1) )
          else
             dens_org(1,i,j) = ( pres_org(2,i,j) + GRAV * dens_org(2,i,j) * cz_org(2,i,j) * 0.5_RP ) &
                             / ( Rdry * temp_org(1) - GRAV * cz_org(2,i,j) * 0.5_RP )
             pres_org(1,i,j) = dens_org(1,i,j) * Rdry * temp_org(1)
             pott_org(1,i,j) = temp_org(1) * ( P00 / pres_org(1,i,j) )**(Rdry/CPdry)
          end if
       else
          dens_org(1:2,i,j) = UNDEF
          pres_org(1:2,i,j) = UNDEF
          pott_org(1:2,i,j) = UNDEF
       end if

    end do
    end do

    return
  end subroutine ParentAtmosInputSCALE

  !-----------------------------------------------------------------------------
  !> Land Setup
  subroutine ParentLandSetupSCALE( &
       ldims )
    use scale_comm_cartesC_nest, only: &
       PARENT_IMAX,                                     &
       PARENT_JMAX,                                     &
       PARENT_LKMAX,                                    &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y => COMM_CARTESC_NEST_TILE_NUM_Y
    implicit none

    integer, intent(out) :: ldims(3)
    !---------------------------------------------------------------------------

    LOG_INFO("ParentLandSetupSCALE",*) 'Real Case/Land Input File Type: SCALE-RM'

    ldims(1) = PARENT_LKMAX(handle)
    ldims(2) = PARENT_IMAX (handle) * NEST_TILE_NUM_X
    ldims(3) = PARENT_JMAX (handle) * NEST_TILE_NUM_Y

    if ( .NOT. allocated(read2D) ) then
       allocate( read2D(PARENT_IMAX(handle),PARENT_JMAX(handle)) )
    endif

    allocate( read3DL(ldims(1),PARENT_IMAX(handle),PARENT_JMAX(handle)) )

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
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R   => CONST_D2R
    use scale_comm_cartesC_nest, only: &
       PARENT_IMAX,                                     &
       PARENT_JMAX,                                     &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_ID    => COMM_CARTESC_NEST_TILE_ID
    use scale_file, only: &
       FILE_open, &
       FILE_read, &
       FILE_get_dataInfo
    use scale_file_CARTESC, only: &
       FILE_CARTESC_read
    implicit none

    real(RP),         intent(out) :: tg_org   (:,:,:)
    real(RP),         intent(out) :: strg_org (:,:,:)
    real(RP),         intent(out) :: lst_org  (:,:)
    real(RP),         intent(out) :: ust_org  (:,:)
    real(RP),         intent(out) :: albg_org (:,:,:,:)
    real(RP),         intent(out) :: topo_org (:,:)
    real(RP),         intent(out) :: lmask_org(:,:)
    real(RP),         intent(out) :: llon_org (:,:)
    real(RP),         intent(out) :: llat_org (:,:)
    real(RP),         intent(out) :: lz_org   (:)
    character(len=*), intent(in)  :: basename_land
    integer,          intent(in)  :: ldims(3)
    logical,          intent(in)  :: use_file_landwater ! use land water data from files
    integer,          intent(in)  :: it

    integer :: rank
    integer :: fid
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye
    logical :: existed

    integer :: i
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

       call FILE_open( BASENAME_land,                 & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       call FILE_CARTESC_read( fid, "LAND_TEMP",  read3DL(:,:,:), step=it )
       tg_org(1:ldims(1),xs:xe,ys:ye) = read3DL(:,:,:)

       if( use_file_landwater )then
          call FILE_CARTESC_read( fid, "LAND_WATER", read3DL(:,:,:), step=it )
          strg_org(1:ldims(1),xs:xe,ys:ye) = read3DL(:,:,:)
       endif

       call FILE_CARTESC_read( fid, "lon", read2D(:,:) )
       llon_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FILE_CARTESC_read( fid, "lat", read2D(:,:) )
       llat_org (xs:xe,ys:ye)  = read2D(:,:) * D2R

       call FILE_CARTESC_read( fid, "LAND_SFC_TEMP", read2D(:,:), step=it )
       lst_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_get_dataInfo( fid, "URBAN_SFC_TEMP", istep=it, existed=existed )
       if ( existed ) then
          call FILE_CARTESC_read( fid, "URBAN_SFC_TEMP", read2D(:,:), step=it )
          ust_org(xs:xe,ys:ye) = read2D(:,:)
       else
          ust_org(xs:xe,ys:ye) = UNDEF
       end if

       call FILE_CARTESC_read( fid, "LAND_SFC_ALB_IR_dir", read2D(:,:), step=it )
       albg_org(xs:xe,ys:ye,I_R_direct ,I_R_IR ) = read2D(:,:)

       call FILE_CARTESC_read( fid, "LAND_SFC_ALB_IR_dif", read2D(:,:), step=it )
       albg_org(xs:xe,ys:ye,I_R_diffuse,I_R_IR ) = read2D(:,:)

       call FILE_CARTESC_read( fid, "LAND_SFC_ALB_NIR_dir", read2D(:,:), step=it )
       albg_org(xs:xe,ys:ye,I_R_direct ,I_R_NIR) = read2D(:,:)

       call FILE_CARTESC_read( fid, "LAND_SFC_ALB_NIR_dif", read2D(:,:), step=it )
       albg_org(xs:xe,ys:ye,I_R_diffuse,I_R_NIR) = read2D(:,:)

       call FILE_CARTESC_read( fid, "LAND_SFC_ALB_VIS_dir", read2D(:,:), step=it )
       albg_org(xs:xe,ys:ye,I_R_direct ,I_R_VIS) = read2D(:,:)

       call FILE_CARTESC_read( fid, "LAND_SFC_ALB_VIS_dif", read2D(:,:), step=it )
       albg_org(xs:xe,ys:ye,I_R_diffuse,I_R_VIS) = read2D(:,:)

       call FILE_CARTESC_read( fid, "topo", read2D(:,:) )
       topo_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_CARTESC_read( fid, "lsmask", read2D(:,:) )
       lmask_org(xs:xe,ys:ye) = read2D(:,:)

    end do

    call FILE_read( fid, "lz", lz_org(:) )

    return
  end subroutine ParentLandInputSCALE

  !-----------------------------------------------------------------------------
  !> Ocean Setup
  subroutine ParentOceanSetupSCALE( &
       odims )
    use scale_comm_cartesC_nest, only: &
       PARENT_IMAX,                                     &
       PARENT_JMAX,                                     &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y => COMM_CARTESC_NEST_TILE_NUM_Y
    implicit none

    integer, intent(out) :: odims(2)
    !---------------------------------------------------------------------------

    LOG_INFO("ParentOceanSetupSCALE",*) 'Real Case/Ocean Input File Type: SCALE-RM'

    odims(1) = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    odims(2) = PARENT_JMAX(handle) * NEST_TILE_NUM_Y

    if ( .NOT. allocated(read2D) ) then
       allocate( read2D(PARENT_IMAX(handle),PARENT_JMAX(handle)) )
    end if

    allocate( read3DO(1,PARENT_IMAX(handle),PARENT_JMAX(handle)) )

    return
  end subroutine ParentOceanSetupSCALE

  !-----------------------------------------------------------------------------
  subroutine ParentOceanOpenSCALE( &
       olon_org,      &
       olat_org,      &
       omask_org,     &
       basename_ocean )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_comm_cartesC_nest, only: &
       PARENT_IMAX,                                     &
       PARENT_JMAX,                                     &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_ID    => COMM_CARTESC_NEST_TILE_ID
    use scale_file, only: &
       FILE_open
    use scale_file_CARTESC, only: &
       FILE_CARTESC_read
    implicit none

    real(RP),         intent(out) :: olon_org (:,:)
    real(RP),         intent(out) :: olat_org (:,:)
    real(RP),         intent(out) :: omask_org(:,:)
    character(len=*), intent(in)  :: basename_ocean

    integer :: rank
    integer :: xloc, yloc
    integer :: xs, xe, ys, ye

    integer :: fid
    integer :: i
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

       call FILE_open( BASENAME_ocean,                & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       call FILE_CARTESC_read( fid, "lon", read2D(:,:) )
       olon_org (xs:xe,ys:ye) = read2D(:,:) * D2R

       call FILE_CARTESC_read( fid, "lat", read2D(:,:) )
       olat_org (xs:xe,ys:ye) = read2D(:,:) * D2R

       call FILE_CARTESC_read( fid, "lsmask", read2D(:,:) )
       omask_org(xs:xe,ys:ye) = read2D(:,:)

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
       it              )
    use scale_comm_cartesC_nest, only: &
       PARENT_IMAX,                                     &
       PARENT_JMAX,                                     &
       NEST_TILE_NUM_X => COMM_CARTESC_NEST_TILE_NUM_X, &
       NEST_TILE_ID    => COMM_CARTESC_NEST_TILE_ID
    use scale_file, only: &
       FILE_open, &
       FILE_get_dataInfo
    use scale_file_CARTESC, only: &
       FILE_CARTESC_read
    implicit none

    real(RP),         intent(out) :: tw_org   (:,:)
    real(RP),         intent(out) :: sst_org  (:,:)
    real(RP),         intent(out) :: albw_org (:,:,:,:)
    real(RP),         intent(out) :: z0w_org  (:,:)
    real(RP),         intent(out) :: omask_org(:,:)
    character(len=*), intent(in)  :: basename_ocean
    integer,          intent(in)  :: it

    integer :: rank

    integer :: fid
    integer :: ndim
    integer :: i
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye
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

       call FILE_open( BASENAME_ocean,                & ! (in)
                       fid,                           & ! (out)
                       aggregate=.false., rankid=rank ) ! (in)

       call FILE_get_dataInfo( fid, "OCEAN_TEMP", dim_rank=ndim )
       select case (ndim)
       case ( 2 ) ! for old file
          call FILE_CARTESC_read( fid, "OCEAN_TEMP", read2D(:,:), step=it )
          tw_org(xs:xe,ys:ye) = read2D(:,:)
       case ( 3 )
          call FILE_CARTESC_read( fid, "OCEAN_TEMP", read3DO(:,:,:), step=it )
          tw_org(xs:xe,ys:ye) = read3DO(1,:,:)
       end select

       call FILE_CARTESC_read( fid, "OCEAN_SFC_TEMP", read2D(:,:), step=it )
       sst_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_CARTESC_read( fid, "OCEAN_SFC_ALB_IR_dir", read2D(:,:), step=it )
       albw_org(xs:xe,ys:ye,I_R_direct ,I_R_IR ) = read2D(:,:)

       call FILE_CARTESC_read( fid, "OCEAN_SFC_ALB_IR_dif", read2D(:,:), step=it )
       albw_org(xs:xe,ys:ye,I_R_diffuse,I_R_IR ) = read2D(:,:)

       call FILE_CARTESC_read( fid, "OCEAN_SFC_ALB_NIR_dir", read2D(:,:), step=it )
       albw_org(xs:xe,ys:ye,I_R_direct ,I_R_NIR) = read2D(:,:)

       call FILE_CARTESC_read( fid, "OCEAN_SFC_ALB_NIR_dif", read2D(:,:), step=it )
       albw_org(xs:xe,ys:ye,I_R_diffuse,I_R_NIR) = read2D(:,:)

       call FILE_CARTESC_read( fid, "OCEAN_SFC_ALB_VIS_dir", read2D(:,:), step=it )
       albw_org(xs:xe,ys:ye,I_R_direct ,I_R_VIS) = read2D(:,:)

       call FILE_CARTESC_read( fid, "OCEAN_SFC_ALB_VIS_dif", read2D(:,:), step=it )
       albw_org(xs:xe,ys:ye,I_R_diffuse,I_R_VIS) = read2D(:,:)

       call FILE_CARTESC_read( fid, "OCEAN_SFC_Z0M", read2D(:,:), step=it )
       z0w_org(xs:xe,ys:ye) = read2D(:,:)

       call FILE_CARTESC_read( fid, "lsmask", read2D(:,:) )
       omask_org(xs:xe,ys:ye) = read2D(:,:)

    end do

    return
  end subroutine ParentOceanInputSCALE

end module mod_realinput_scale
