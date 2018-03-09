!-------------------------------------------------------------------------------
!> module REAL input NICAM
!!
!! @par Description
!!          read data from NICAM file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-05-24 (S.Nishizawa)   [new] split from mod_realinput.f90
!!
!<
module mod_realinput_nicam
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

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtmosSetupNICAM
  public :: ParentAtmosOpenNICAM
  public :: ParentAtmosInputNICAM
  public :: ParentLandSetupNICAM
  public :: ParentLandInputNICAM
  public :: ParentOceanSetupNICAM
  public :: ParentOceanOpenNICAM
  public :: ParentOceanInputNICAM

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
  real(RP), parameter :: missval_tg   = 50.0_RP
  real(RP), parameter :: missval_strg = 0.0_RP
  real(RP), parameter :: missval_sst  = 273.1506_RP


  real(SP), allocatable :: read1DX(:)
  real(SP), allocatable :: read1DY(:)
  real(SP), allocatable :: read1DZ(:)
  real(SP), allocatable :: read3D (:,:,:)

  real(SP), allocatable :: read3DL(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtmosSetupNICAM( &
      dims,    &
      timelen, &
      basename_org )
    use scale_file, only: &
         FILE_Get_Shape
    implicit none

    integer,          intent(out) :: dims(6)
    integer,          intent(out) :: timelen
    character(len=*), intent(in)  :: basename_org

    character(len=H_LONG) :: basename
    integer :: dims_ncm(4)

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atmos Input File Type: NICAM-NETCDF'
    basename = "ms_pres"//trim(basename_org)
    call FILE_Get_Shape( basename, "ms_pres", & ! (in)
                         dims_ncm(:),         & ! (out)
                         single=.true.        ) ! (in)
    timelen = dims_ncm(4)

    ! full level
    dims(1) = dims_ncm(3)
    dims(2) = dims_ncm(1)
    dims(3) = dims_ncm(2)
    ! half level
    ! nicam lat-lon data doesn't have staggered grid system
    dims(4) = dims(1)
    dims(5) = dims(2)
    dims(6) = dims(3)

    allocate( read1DX( dims(2)                   ) )
    allocate( read1DY( dims(3)                   ) )
    allocate( read1DZ( dims(1)                   ) )
    allocate( read3D ( dims(1), dims(2), dims(3) ) )

    return
  end subroutine ParentAtmosSetupNICAM

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosOpenNICAM( &
       lon_org,      &
       lat_org,      &
       cz_org,       &
       basename_num, &
       dims          )
    use scale_const, only: &
         D2R => CONST_D2R
    use scale_file, only: &
         FILE_open, &
         FILE_Read
    implicit none

    real(RP),         intent(out) :: lon_org(:,:)
    real(RP),         intent(out) :: lat_org(:,:)
    real(RP),         intent(out) :: cz_org(:,:,:)
    character(len=*), intent(in)  :: basename_num
    integer,          intent(in)  :: dims(6)

    character(len=H_LONG) :: basename

    integer :: fid
    integer :: k, i, j

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[AtmosOpenNICAM]'

    basename = "ms_pres"//trim(basename_num)
    call FILE_open( basename, fid, single=.true., postfix="" )

    call FILE_read( fid, "lon", read1DX(:) )
    do j = 1, dims(3)
       lon_org (:,j)  = read1DX(:) * D2R
    enddo

    call FILE_read( fid, "lat", read1DY(:) )
    do i = 1, dims(2)
       lat_org (i,:)  = read1DY(:) * D2R
    enddo

    call FILE_read( fid, "lev", read1DZ(:) )
    do j = 1, dims(3)
    do i = 1, dims(2)
       cz_org(3:,i,j) = read1DZ(:)
    end do
    end do

    call FILE_read( fid, "ms_pres", read3D(:,:,:) )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = dims(1)+2, 3, -1
          ! missing data implies under ground
          ! For more accurate surface height, we need topograph data
          ! So surface data is not used at this moment
          if ( read3D(i,j,k-2) <= 0.0_RP ) then
!             cz_org(k,i,j) = max( ( cz_org(k+1,i,j) + cz_org(k,i,j) ) * 0.5_RP, 0.0_RP )
             cz_org(1:k,i,j) = 0.0_RP
             exit
          else
             cz_org(k,i,j) = read3D(i,j,k)
          endif
       enddo
    enddo
    enddo

    return
  end subroutine ParentAtmosOpenNICAM

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosInputNICAM( &
       velz_org,     &
       velx_org,     &
       vely_org,     &
       pres_org,     &
       temp_org,     &
       qtrc_org,     &
       basename_num, &
       dims,         &
       it            )
    use scale_const, only: &
       Rdry => CONST_Rdry, &
       CPdry => CONST_CPdry, &
       P00 => CONST_PRE00, &
       D2R => CONST_D2R, &
       EPS => CONST_EPS
    use scale_file, only: &
       FILE_read
    use scale_atmos_hydrometeor, only: &
       I_QV
    implicit none

    real(RP),         intent(out) :: velz_org(:,:,:)
    real(RP),         intent(out) :: velx_org(:,:,:)
    real(RP),         intent(out) :: vely_org(:,:,:)
    real(RP),         intent(out) :: pres_org(:,:,:)
    real(RP),         intent(out) :: temp_org(:,:,:)
    real(RP),         intent(out) :: qtrc_org(:,:,:,:)
    character(len=*), intent(in)  :: basename_num
    integer,          intent(in)  :: dims(6)
    integer,          intent(in)  :: it

    real(RP) :: tsfc_org (dims(2),dims(3))
    real(RP) :: slp_org  (dims(2),dims(3))
    real(RP) :: qvsfc_org(dims(2),dims(3))

    real(RP) :: pott
    real(RP) :: RovCP
    real(RP) :: CPovR

    integer :: k, i, j

    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[AtmosInputNICAM]'

    basename = "ms_u"//trim(basename_num)
    call FILE_read( basename, "ms_u", read3D(:,:,:), step=it, rankid=myrank, single=.true., postfix="" )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 1, dims(1)
          if ( abs(abs(read3D(i,j,k))-300.0_SP) < sqrt(EPS) ) then ! missing value
             velx_org(k+2,i,j) = 0.0_RP
          else
             velx_org(k+2,i,j) = read3D(i,j,k)
          end if
       end do
       velx_org(1:2,i,j) = 0.0_RP
    end do
    end do

    basename = "ms_v"//trim(basename_num)
    call FILE_read( basename, "ms_v", read3D(:,:,:), step=it, rankid=myrank, single=.true., postfix="" )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 1, dims(1)
          if ( abs(abs(read3D(i,j,k))-300.0_SP) < sqrt(EPS) ) then ! missing value
             vely_org(k+2,i,j) = 0.0_RP
          else
             vely_org(k+2,i,j) = read3D(i,j,k)
          end if
       end do
       vely_org(1:2,i,j) = 0.0_RP
    end do
    end do

    velz_org(:,:,:) = 0.0_RP !> cold initialize for vertical velocity


    ! The surface height is not available, so t2 is not used, at this moment.
!    basename = "ss_t2m"//trim(basename_num)
!    call FILE_read( basename, "ss_t2m", tsfc_org(:,:), step=it, rankid=myrank, single=.true., postfix="" )

    basename = "ms_tem"//trim(basename_num)
    call FILE_read( basename, "ms_tem", read3D(:,:,:), step=it, rankid=myrank, single=.true., postfix="" )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = dims(1), 1, -1
          if ( read3D(i,j,k) <= 50.0_RP ) then ! missing value
             temp_org(1:k+2,i,j) = temp_org(k+3,i,j)
             exit
          else
             temp_org(k+2,i,j) = read3D(i,j,k)
          end if
!          if( k==1 ) temp_org(2,i,j) = tsfc_org(i,j) ! no missing value case
       end do
    end do
    end do

    RovCP = Rdry/CPdry
    CPovR = CPdry/Rdry

    basename = "ss_slp"//trim(basename_num)
    call FILE_read( basename, "ss_slp", slp_org(:,:), step=it, rankid=myrank, single=.true., postfix="" )

    basename = "ms_pres"//trim(basename_num)
    call FILE_read( basename, "ms_pres", read3D(:,:,:), step=it, rankid=myrank, single=.true., postfix="" )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = dims(1), 1, -1
          if ( read3D(i,j,k) < 0.0_RP ) then ! missing value
             exit
          else
             pres_org(k+2,i,j) = read3D(i,j,k)
          end if
       end do
       ! If data has missing value, k is the top layer having missing value,
       ! otherwise k is zero.
       pott = temp_org(k+3,i,j) * (P00/pres_org(k+3,i,j))**RovCP ! lowest level
!       pres_org(k+2,i,j) = P00 * (temp_org(k+2,i,j)/pott)**CPovR ! surface
       pres_org(1:k+2,i,j) = slp_org(i,j)                        ! sea level
       temp_org(1:k+2,i,j) = pott * (slp_org(i,j)/P00)**RovCP    ! sea level
    end do
    end do

    basename = "ss_q2m"//trim(basename_num)
    call FILE_read( basename, "ss_q2m", qvsfc_org(:,:), step=it, rankid=myrank, single=.true., postfix="" )

    qtrc_org(:,:,:,:) = 0.0_RP
    basename = "ms_qv"//trim(basename_num)
    call FILE_read( basename, "ms_qv", read3D(:,:,:), step=it, rankid=myrank, single=.true., postfix="" )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = dims(1), 1, -1
          if ( read3D(i,j,k) < 0.0_RP ) then ! missing value
             exit
          else
             qtrc_org(k+2,i,j,I_QV) = read3D(i,j,k)
          end if
       end do
       qtrc_org(1:k+2,i,j,I_QV) = qvsfc_org(i,j)
    end do
    end do


    return
  end subroutine ParentAtmosInputNICAM

  !-----------------------------------------------------------------------------
  !> Land Setup
  subroutine ParentLandSetupNICAM( &
       ldims,    &
       basename_org )
    use scale_file, only: &
         FILE_Get_Shape
    implicit none

    integer,          intent(out) :: ldims(3)
    character(len=*), intent(in)  :: basename_org

    character(len=H_LONG) :: basename
    integer :: dims_ncm(4)

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Land Input File Type: NICAM-NETCDF'
    basename = "la_tg"//trim(basename_org)
    call FILE_Get_Shape( basename, "la_tg", & ! (in)
                         dims_ncm(:),       & ! (out)
                         single=.true.      ) ! (in)
    ! land
    ldims(1) = dims_ncm(3) ! vertical grid of land model
    ldims(2) = dims_ncm(1)
    ldims(3) = dims_ncm(2)

    if ( .not. allocated(read1DX) ) then
       allocate( read1DX( ldims(2)                       ) )
       allocate( read1DY( ldims(3)                       ) )
    end if

    allocate( read3DL ( ldims(2), ldims(3), ldims(1) ) )

    return
  end subroutine ParentLandSetupNICAM

  !-----------------------------------------------------------------------------
  subroutine ParentLandInputNICAM( &
      tg_org,             &
      strg_org,           &
      lst_org,            &
      llon_org,           &
      llat_org,           &
      lz_org,             &
      topo_org,           &
      lmask_org,          &
      basename_num,       &
      ldims,              &
      use_file_landwater, &
      it                  )
    use scale_const, only: &
         UNDEF => CONST_UNDEF, &
         D2R   => CONST_D2R,   &
         TEM00 => CONST_TEM00, &
         EPS   => CONST_EPS
    use scale_file, only: &
         FILE_Read
    implicit none

    real(RP), intent(out) :: tg_org(:,:,:)
    real(RP), intent(out) :: strg_org(:,:,:)
    real(RP), intent(out) :: lst_org(:,:)
    real(RP), intent(out) :: llon_org(:,:)
    real(RP), intent(out) :: llat_org(:,:)
    real(RP), intent(out) :: lz_org(:)
    real(RP), intent(out) :: topo_org(:,:)
    real(RP), intent(out) :: lmask_org(:,:)
    character(len=*), intent(in) :: basename_num
    integer,          intent(in) :: ldims(3)
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    integer,          intent(in) :: it

    ! work
    integer :: k, i, j

    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[LandInputNICAM]'

    basename = "la_tg"//trim(basename_num)
    call FILE_Read( basename, "lev", lz_org(:), single=.true., postfix="" )

    call FILE_Read( basename, "lon", read1DX(:), single=.true., postfix="" )
    do j = 1, ldims(3)
       llon_org (:,j)  = read1DX(:) * D2R
    enddo

    call FILE_Read( basename, "lat", read1DY(:), single=.true., postfix="" )
    do i = 1, ldims(2)
       llat_org (i,:)  = read1DY(:) * D2R
    enddo

    basename = "lsmask"//trim(basename_num)
    call FILE_read( basename, "lsmask", lmask_org(:,:), step=it, rankid=myrank, single=.true., allow_missing=.true., missing_value=UNDEF, postfix="" )

    basename = "la_tg"//trim(basename_num)
    call FILE_read( basename, "la_tg", read3DL(:,:,:), step=it, rankid=myrank, single=.true., postfix="" )
    do j = 1, ldims(3)
    do i = 1, ldims(2)
    do k = 1, ldims(1)
       tg_org(k,i,j) = read3DL(i,j,k)
    end do
    end do
    end do

    if( use_file_landwater ) then
       basename = "la_wg"//trim(basename_num)
       call FILE_read( basename, "la_wg", read3DL(:,:,:), step=it, rankid=myrank, single=.true., postfix="" )
       do j = 1, ldims(3)
       do i = 1, ldims(2)
       do k = 1, ldims(1)
          strg_org(k,i,j) = read3DL(i,j,k)
       end do
       end do
       end do
    end if

    basename = "ss_tem_sfc"//trim(basename_num)
    call FILE_read( basename, "ss_tem_sfc", lst_org(:,:), step=it, rankid=myrank, single=.true., postfix="" )

    ! replace missing value
    do j = 1, ldims(3)
    do i = 1, ldims(2)
    do k = 1, ldims(1)
       if ( abs(tg_org(k,i,j)  -missval_tg  ) < EPS ) tg_org(k,i,j)   = UNDEF
       if ( abs(strg_org(k,i,j)-missval_strg) < EPS ) strg_org(k,i,j) = UNDEF
    end do
    end do
    end do

    topo_org = UNDEF

    return
  end subroutine ParentLandInputNICAM

  !-----------------------------------------------------------------------------
  !> Ocean Setup
  subroutine ParentOceanSetupNICAM( &
       odims,    &
       timelen, &
       basename_org )
    use scale_file, only: &
         FILE_Get_Shape
    implicit none

    integer,          intent(out) :: odims(2)
    integer,          intent(out) :: timelen
    character(len=*), intent(in)  :: basename_org

    character(len=H_LONG) :: basename
    integer :: dims_ncm(4)

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Ocean Input File Type: NICAM-NETCDF'

    basename = "oa_sst"//trim(basename_org)
    call FILE_Get_Shape( basename, "oa_sst", &
                         dims_ncm(:),        &
                         single=.true.       )
    odims(1) = dims_ncm(1)
    odims(2) = dims_ncm(2)

    timelen = dims_ncm(4)

    if ( .not. allocated(read1DX) ) then
       allocate( read1DX( odims(1)                       ) )
       allocate( read1DY( odims(2)                       ) )
    end if

    return
  end subroutine ParentOceanSetupNICAM

  !-----------------------------------------------------------------------------
  subroutine ParentOceanOpenNICAM( &
       olon_org,     &
       olat_org,     &
       omask_org,    &
       basename_num, &
       odims         )
    use scale_const, only: &
         D2R   => CONST_D2R, &
         UNDEF => CONST_UNDEF
    use scale_file, only: &
         FILE_read
    implicit none

    real(RP),         intent(out) :: olon_org(:,:)
    real(RP),         intent(out) :: olat_org(:,:)
    real(RP),         intent(out) :: omask_org(:,:)
    character(len=*), intent(in)  :: basename_num
    integer,          intent(in)  :: odims(2)

    character(len=H_LONG) :: basename

    integer :: k, i, j

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OceanOpenNICAM]'

    basename = "oa_sst"//trim(basename_num)
    call FILE_read( basename, "lon", read1DX(:), single=.true., postfix="" )
    do j = 1, odims(2)
       olon_org (:,j)  = read1DX(:) * D2R
    enddo

    call FILE_read( basename, "lat", read1DY(:), single=.true., postfix="" )
    do i = 1, odims(1)
       olat_org (i,:)  = read1DY(:) * D2R
    enddo

    basename = "lsmask"//trim(basename_num)
    call FILE_read( basename, "lsmask", omask_org(:,:), step=1, rankid=myrank, single=.true., allow_missing=.true., missing_value=UNDEF, postfix="" )

    return
  end subroutine ParentOceanOpenNICAM

  !-----------------------------------------------------------------------------
  subroutine ParentOceanInputNICAM( &
       tw_org,       &
       sst_org,      &
       basename_num, &
       odims,        &
       omask_org,    &
       it            )
    use scale_const, only: &
         UNDEF => CONST_UNDEF, &
         D2R   => CONST_D2R,   &
         TEM00 => CONST_TEM00, &
         EPS   => CONST_EPS
    use scale_file, only: &
         FILE_read
    implicit none

    real(RP), intent(out) :: tw_org(:,:)
    real(RP), intent(out) :: sst_org(:,:)
    character(len=*), intent(in) :: basename_num
    integer,          intent(in) :: odims(2)
    real(RP),         intent(in) :: omask_org(:,:)
    integer,          intent(in) :: it

    ! work
    integer :: i, j

    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OceanInputNICAM]'


    basename = "oa_sst"//trim(basename_num)
    call FILE_read( basename, "oa_sst", sst_org(:,:), step=1, rankid=myrank, single=.true., postfix="" )

    !basename = "oa_ice"//trim(basename_num)
    !call FILE_read( basename, "oa_ice", ice_org(:,:), step=1, rankid=myrank, single=.true., postfix="" )


    ! SST: retrieve SST data around coast
    do j = 1, odims(2)
    do i = 1, odims(1)
       if ( abs(omask_org(i,j)-1.0_RP) < EPS ) then ! land data
          cycle
       else
          sst_org(i,j) = ( sst_org(i,j) - missval_sst*omask_org(i,j) ) &
                       / ( 1.0_RP - omask_org(i,j) )
       end if
    end do
    end do

    tw_org = sst_org

    return
  end subroutine ParentOceanInputNICAM

end module mod_realinput_nicam
