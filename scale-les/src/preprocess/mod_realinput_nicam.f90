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
     PRC_master,            &
     PRC_MPIstop
  use scale_external_io, only: &
     iNICAM

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtomSetupNICAM
  public :: ParentAtomOpenNICAM
  public :: ParentAtomInputNICAM
  public :: ParentSurfaceInputNICAM

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
  real(SP), allocatable :: read3DS(:,:,:,:)
  real(SP), allocatable :: read4D (:,:,:,:)

  real(SP), allocatable :: read1DLZ(:)
  real(SP), allocatable :: read4DL(:,:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ParentAtomSetupNICAM( &
      dims,    &
      timelen, &
      basename_org )
    use gtool_file, only: &
         FileGetShape
    implicit none

    integer,               intent(out) :: dims(11)
    integer,               intent(out) :: timelen
    character(len=H_LONG), intent(in)  :: basename_org

    character(len=H_LONG) :: basename
    integer :: dims_ncm(4)

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: NICAM-NETCDF'
    basename = "ms_pres"//trim(basename_org)
    call FileGetShape( dims_ncm(:), trim(basename), "ms_pres", 1, single=.true. )
    timelen = dims_ncm(4)

    ! full level
    dims(1) = dims_ncm(3)
    dims(2) = dims_ncm(1)
    dims(3) = dims_ncm(2)
    ! half level
    dims(4) = dims_ncm(1) ! nicam lat-lon data doesn't have staggered grid system
    dims(5) = dims_ncm(2) ! nicam lat-lon data doesn't have staggered grid system
    dims(6) = dims_ncm(3) ! nicam lat-lon data doesn't have staggered grid system
    basename = "la_tg"//trim(basename_org)
    call FileGetShape( dims_ncm(:), trim(basename), "la_tg", 1, single=.true. )
    ! land
    dims(7) = dims_ncm(3) ! vertical grid of land model [tentative]
    dims(8) = dims(2)
    dims(9) = dims(3)
    ! sst
    dims(10) = dims(2)
    dims(11) = dims(3)

    allocate( read1DX( dims(2)                      ) )
    allocate( read1DY( dims(3)                      ) )
    allocate( read1DZ( dims(1)                      ) )
    allocate( read3DS( 1,       dims(2), dims(3), 1 ) )
    allocate( read4D ( dims(1), dims(2), dims(3), 1 ) )

    allocate( read1DLZ( dims(7)                      ) )
    allocate( read4DL ( dims(7), dims(8), dims(9), 1 ) )


    return
  end subroutine ParentAtomSetupNICAM

  !-----------------------------------------------------------------------------
  subroutine ParentAtomOpenNICAM( &
       lon_org,      &
       lat_org,      &
       cz_org,       &
       basename_num, &
       dims          )
    use scale_const, only: &
         D2R => CONST_D2R
    use gtool_file, only: &
         FileRead
    use scale_external_io, only: &
         ExternalFileRead
    implicit none

    real(RP), intent(out) :: lon_org(:,:)
    real(RP), intent(out) :: lat_org(:,:)
    real(RP), intent(out) :: cz_org(:,:,:)
    character(len=H_LONG), intent(in) :: basename_num
    integer,  intent(in)  :: dims(9)

    character(len=H_LONG) :: basename

    integer :: k, i, j

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OpenNICAM]'

    basename = "ms_pres"//trim(basename_num)
    call FileRead( read1DX(:), trim(basename), "lon", 1, 1, single=.true. )
    do j = 1, dims(2)
       lon_org (:,j)  = read1DX(:) * D2R
    enddo

    call FileRead( read1DY(:), trim(basename), "lat", 1, 1, single=.true. )
    do i = 1, dims(1)
       lat_org (i,:)  = read1DY(:) * D2R
    enddo

    call FileRead( read1DZ(:), trim(basename), "lev", 1, 1, single=.true. )
    do j = 1, dims(3)
    do i = 1, dims(2)
       cz_org(3:,i,j) = read1DZ(:)
    end do
    end do

    call ExternalFileRead( read4D(:,:,:,:), trim(basename), "ms_pres", 1, 1, myrank, iNICAM, single=.true. )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = dims(1)+2, 3, -1
          ! missing data implies under ground
          ! For more accurate surface height, we need topograph data
          if ( read4D(k-2,i,j,1) <= 0.0_RP ) then
             cz_org(k,i,j) = ( cz_org(k+1,i,j) + cz_org(k,i,j) ) * 0.5_RP
             cz_org(k-1,i,j) = 0.0_RP
             cz_org(1:k-2,i,j) = -1e5_RP
             exit
          endif
       enddo
    enddo
    enddo

    return
  end subroutine ParentAtomOpenNICAM
  !-----------------------------------------------------------------------------
  subroutine ParentAtomInputNICAM( &
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
    use scale_external_io, only: &
         ExternalFileRead, &
         ExternalFileReadOffset
    implicit none

    real(RP),         intent(out) :: velz_org(:,:,:)
    real(RP),         intent(out) :: velx_org(:,:,:)
    real(RP),         intent(out) :: vely_org(:,:,:)
    real(RP),         intent(out) :: pres_org(:,:,:)
    real(RP),         intent(out) :: temp_org(:,:,:)
    real(RP),         intent(out) :: qtrc_org(:,:,:,:)
    character(LEN=*), intent(in)  :: basename_num
    integer,          intent(in)  :: dims(7)
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

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputNICAM]'

    !> [scale-offset]
    basename = "ms_u"//trim(basename_num)
    call ExternalFileReadOffset( read4D(:,:,:,:), trim(basename), "ms_u", it, it, myrank, iNICAM, single=.true. )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 1, dims(1)
          if ( abs(abs(read4D(k,i,j,1))-300.0_SP) < sqrt(EPS) ) then ! missing value
             velx_org(k+2,i,j) = 0.0_RP
          else
             velx_org(k+2,i,j) = real( read4D(k,i,j,1), kind=RP )
          end if
       end do
       velx_org(1:2,i,j) = 0.0_RP
    end do
    end do

    basename = "ms_v"//trim(basename_num)
    call ExternalFileReadOffset( read4D(:,:,:,:),  &
                                 trim(basename),   &
                                 "ms_v",           &
                                 it, it,      &
                                 myrank,           &
                                 iNICAM,           &
                                 single=.true.     )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 1, dims(1)
          if ( abs(abs(read4D(k,i,j,1))-300.0_SP) < sqrt(EPS) ) then ! missing value
             vely_org(k+2,i,j) = 0.0_RP
          else
             vely_org(k+2,i,j) = real( read4D(k,i,j,1), kind=RP )
          end if
       end do
       vely_org(1:2,i,j) = 0.0_RP
    end do
    end do

    velz_org(:,:,:) = 0.0_RP !> cold initialize for vertical velocity


    basename = "ss_t2m"//trim(basename_num)
    call ExternalFileRead( read3DS(:,:,:,:), trim(basename), "ss_t2m", it, it, myrank, iNICAM, single=.true. )
    tsfc_org(:,:) = real( read3DS(1,:,:,1), kind=RP )
    basename = "ms_tem"//trim(basename_num)
    call ExternalFileReadOffset( read4D(:,:,:,:), trim(basename), "ms_tem", it, it, myrank, iNICAM, single=.true. )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = dims(1), 1, -1
          if ( read4D(k,i,j,1) <= 50.0_SP ) then ! missing value
             temp_org(k+2,i,j) = tsfc_org(i,j)
             exit
          else
             temp_org(k+2,i,j) = real( read4D(k,i,j,1), kind=RP )
          end if
          if( k==-1 ) temp_org(2,i,j) = tsfc_org(i,j) ! no missing value case
       end do
    end do
    end do

    RovCP = Rdry/CPdry
    CPovR = CPdry/Rdry

    basename = "ss_slp"//trim(basename_num)
    call ExternalFileReadOffset( read3DS(:,:,:,:), trim(basename), "ss_slp", it, it, myrank, iNICAM, single=.true. )
    slp_org(:,:) = real( read3DS(1,:,:,1), kind=RP )
    basename = "ms_pres"//trim(basename_num)
    call ExternalFileRead( read4D(:,:,:,:), trim(basename), "ms_pres", it, it, myrank, iNICAM, single=.true. )
    pres_org(3:,:,:) = read4D(:,:,:,1)
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = dims(1), 1, -1
          if ( read4D(k,i,j,1) < 0.0_SP ) then ! missing value
             exit
          else
             pres_org(k+2,i,j) = real( read4D(k,i,j,1), kind=RP )
          end if
       end do
       ! If data has missing value, k is the top layer having missing value,
       ! otherwise k is zero.
       pott = temp_org(k+3,i,j) * (P00/pres_org(k+3,i,j))**RovCP ! lowest level
       pres_org(k+2,i,j) = P00 * (temp_org(k+2,i,j)/pott)**CPovR ! surface
       temp_org(k+1,i,j) = pott * (pres_org(k+1,i,j)/P00)**RovCP ! sea level
       pres_org(k+1,i,j) = slp_org(i,j)                          ! sea level
    end do
    end do


    basename = "ss_q2m"//trim(basename_num)
    call ExternalFileRead( read3DS(:,:,:,:), trim(basename), "ss_q2m", it, it, myrank, iNICAM, single=.true. )
    qvsfc_org(:,:) = real( read3DS(1,:,:,1), kind=RP )
    qtrc_org(:,:,:,:) = 0.0_RP
    basename = "ms_qv"//trim(basename_num)
    call ExternalFileRead( read4D(:,:,:,:), trim(basename), "ms_qv", it, it, myrank, iNICAM, single=.true. )
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = dims(1), 1, -1
          if ( read4D(k,i,j,1) < 0.0_SP ) then ! missing value
             qtrc_org(k+1:k+2,i,j,I_QV) = qvsfc_org(i,j) ! surface and sealevel
             exit
          else
             qtrc_org(k+2,i,j,I_QV) = real( read4D(k,i,j,1), kind=RP )
          end if
       end do
       qtrc_org(1:2,i,j,I_QV) = qvsfc_org(i,j)
    end do
    end do


    return
  end subroutine ParentAtomInputNICAM

  !-----------------------------------------------------------------------------
  subroutine ParentSurfaceInputNICAM( &
      tg_org,             &
      strg_org,           &
      lst_org,            &
      sst_org,            &
      lsmask_org,         &
      lz_org,             &
      basename_num,       &
      dims,               &
      use_file_landwater, &
      it                  )
    use scale_const, only: &
         UNDEF => CONST_UNDEF, &
         D2R   => CONST_D2R,   &
         TEM00 => CONST_TEM00, &
         EPS   => CONST_EPS
    use gtool_file, only: &
         FileRead
    use scale_external_io, only: &
         ExternalFileRead, & 
         ExternalFileReadOffset
    implicit none

    real(RP), intent(out) :: tg_org(:,:,:)
    real(RP), intent(out) :: strg_org(:,:,:)
    real(RP), intent(out) :: lst_org(:,:)
    real(RP), intent(out) :: sst_org(:,:)
    real(RP), intent(out) :: lsmask_org(:,:)
    real(RP), intent(out) :: lz_org(:)
    character(LEN=*), intent(in) :: basename_num
    integer,          intent(in) :: dims(11)
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    integer,          intent(in) :: it

    ! work
    integer :: k, i, j

    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputNICAM-Surface]'

    basename = "la_tg"//trim(basename_num)
    call FileRead( read1DLZ(:), trim(basename), "lev", 1, 1, single=.true. )
    lz_org(:) = read1DLZ(:)

    basename = "lsmask"//trim(basename_num)
    call ExternalFileRead( read3DS(:,:,:,:), &
                           trim(basename),   &
                           "lsmask",         &
                           it, it, &
                           myrank,           &
                           iNICAM,           &
                           single=.true.     )
    lsmask_org(:,:) = real( read3DS(1,:,:,1), kind=RP )

    basename = "la_tg"//trim(basename_num)
    call ExternalFileReadOffset( read4DL(:,:,:,:), &
                                 trim(basename),   &
                                 "la_tg",          &
                                 it, it,      &
                                 myrank,           &
                                 iNICAM,           &
                                 single=.true.     )
    tg_org(:,:,:) = real( read4DL(:,:,:,1), kind=RP )

    if( use_file_landwater ) then
       basename = "la_wg"//trim(basename_num)
       call ExternalFileReadOffset( read4DL(:,:,:,:), &
                                    trim(basename),   &
                                    "la_wg",          &
                                    it, it,      &
                                    myrank,           &
                                    iNICAM,           &
                                    single=.true.     )
       strg_org(:,:,:) = real( read4DL(:,:,:,1), kind=RP )
    end if

    basename = "ss_tem_sfc"//trim(basename_num)
    call ExternalFileRead( read3DS(:,:,:,:), &
                           trim(basename),   &
                           "ss_tem_sfc",     &
                           it, it,       &
                           myrank,           &
                           iNICAM,           &
                           single=.true.     )
    lst_org(:,:) = real( read3DS(1,:,:,1), kind=RP )

    basename = "oa_sst"//trim(basename_num)
    call ExternalFileReadOffset( read3DS(:,:,:,:), &
                                 trim(basename),   &
                                 "oa_sst",         &
                                 1, 1,      &
                                 myrank,           &
                                 iNICAM,           &
                                 single=.true.     )
    sst_org(:,:) = real( read3DS(1,:,:,1), kind=RP )

    !basename = "oa_ice"//trim(basename_num)
    !call ExternalFileRead( read3DS(:,:,:,:), trim(basename), &
    !                       "oa_ice", 1, 1, myrank, iNICAM, single=.true. )
    !ice_org(:,:) = real( read3DS(1,:,:,1), kind=RP )


    ! replace missing value
    do j = 1, dims(9)
    do i = 1, dims(8)
    do k = 1, dims(7)
       if ( abs(tg_org(k,i,j)  -missval_tg  ) < EPS ) tg_org(k,i,j)   = UNDEF
       if ( abs(strg_org(k,i,j)-missval_strg) < EPS ) strg_org(k,i,j) = UNDEF
    end do
    end do
    end do

    ! SST: retrieve SST data around coast
    do j = 1, dims(11)
    do i = 1, dims(10)
       if ( abs(lsmask_org(i,j)-1.0_RP) < EPS ) then ! land data
          cycle
       else
          sst_org(i,j) = ( sst_org(i,j) - missval_sst*lsmask_org(i,j) ) &
                       / ( 1.0_RP - lsmask_org(i,j) )
       end if
    end do
    end do

    return
  end subroutine ParentSurfaceInputNICAM

end module mod_realinput_nicam
!-------------------------------------------------------------------------------
