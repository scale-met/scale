!-------------------------------------------------------------------------------
!> module REAL input WRF-ARW
!!
!! @par Description
!!          read data from WRF-ARW file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-05-24 (S.Nishizawa)   [new] split from mod_realinput.f90
!!
!<
module mod_realinput_wrfarw
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
  public :: ParentAtmosSetupWRFARW
  public :: ParentAtmosOpenWRFARW
  public :: ParentAtmosInputWRFARW
  public :: ParentLandSetupWRFARW
  public :: ParentLandInputWRFARW
  public :: ParentOceanSetupWRFARW
  public :: ParentOceanOpenWRFARW
  public :: ParentOceanInputWRFARW

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
  ! Defined parameters in WRF

  real(RP), parameter :: t0  = 300.0_RP
  real(RP), parameter :: p0  = 1000.0E+2_RP
  real(RP), parameter :: Rd  = 287.04_RP
  real(RP), parameter :: Cp  = 7.0_RP * Rd / 2.0_RP
  real(RP), parameter :: RCP = Rd / Cp

  integer, parameter :: cosin = 1
  integer, parameter :: sine  = 2

  real(RP), allocatable :: read_xy (:,:)
  real(RP), allocatable :: read_xyz(:,:,:)
  real(RP), allocatable :: read_xyw(:,:,:)
  real(RP), allocatable :: read_xyl(:,:,:)

  real(RP), allocatable :: p_org   (:,:,:)
  real(RP), allocatable :: pb_org  (:,:,:)
  real(RP), allocatable :: ph_org  (:,:,:)
  real(RP), allocatable :: phb_org (:,:,:)

  logical, private :: wrfout = .false. ! file type switch (wrfout or wrfrst)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtmosSetupWRFARW( &
      dims,    &
      timelen, &
      basename_org )
    use scale_file, only: &
       FILE_open, &
       FILE_get_dimLength
    implicit none

    integer,          intent(out) :: dims(6)
    integer,          intent(out) :: timelen
    character(len=*), intent(in)  :: basename_org

    logical :: WRF_FILE_TYPE = .false.   ! wrf filetype: T=wrfout, F=wrfrst

    NAMELIST / PARAM_MKINIT_REAL_WRFARW / &
         WRF_FILE_TYPE

    integer :: fid
    integer :: ierr
    logical :: error
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atmos Input File Type: WRF-ARW'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_WRFARW,iostat=ierr)
    if( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_WRFARW. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_MKINIT_REAL_WRFARW)

    call FILE_open( basename_org, fid, rankid=myrank, single=.true., postfix="" )

    call FILE_get_dimLength( fid, "bottom_top",       dims(1) )
    call FILE_get_dimLength( fid, "west_east",        dims(2) )
    call FILE_get_dimLength( fid, "south_north",      dims(3) )
    call FILE_get_dimLength( fid, "bottom_top_stag",  dims(4) )
    call FILE_get_dimLength( fid, "west_east_stag",   dims(5) )
    call FILE_get_dimLength( fid, "south_north_stag", dims(6) )

    call FILE_get_dimLength( fid, "Time", timelen, error=error )
    if ( error ) call FILE_get_dimLength( fid, "time", timelen, error=error)
    if ( error ) timelen = 0

    if ( wrf_file_type ) then
       wrfout = .true.
       if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF History Output'
    else
       wrfout = .false.
       if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF Restart'
    endif


    allocate( read_xy  (dims(2),dims(3)) )
    allocate( read_xyz (dims(2),dims(3),dims(1)) )
    allocate( read_xyw (dims(2),dims(3),dims(4)) )

    allocate( p_org    (dims(1),dims(2),dims(3)) )
    allocate( pb_org   (dims(1),dims(2),dims(3)) )
    allocate( ph_org   (dims(4),dims(2),dims(3)) )
    allocate( phb_org  (dims(4),dims(2),dims(3)) )

    return
  end subroutine ParentAtmosSetupWRFARW

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosOpenWRFARW
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[AtmosOpenWRFARW]'
    return
  end subroutine ParentAtmosOpenWRFARW

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosInputWRFARW( &
       velz_org,      &
       llvelx_org,    &
       llvely_org,    &
       pres_org,      &
       temp_org,      &
       qtrc_org,      &
       lon_org,       &
       lat_org,       &
       cz_org,        &
       basename,      &
       mptype_parent, &
       dims,          &
       it             ) ! (in)
    use scale_const, only: &
       D2R => CONST_D2R, &
       LAPS => CONST_LAPS, &
       Rdry => CONST_Rdry, &
       GRAV => CONST_GRAV
    use scale_file, only: &
       FILE_open, &
       FILE_read
    use scale_atmos_thermodyn, only: &
       THERMODYN_pott => ATMOS_THERMODYN_pott
    use scale_atmos_hydrometeor, only: &
       QHS, &
       QHE, &
       I_QV, &
       I_QC, &
       I_QR, &
       I_QI, &
       I_QS, &
       I_QG, &
       I_NC, &
       I_NR, &
       I_NI, &
       I_NS, &
       I_NG
    implicit none

    real(RP),         intent(out) :: velz_org(:,:,:)
    real(RP),         intent(out) :: llvelx_org(:,:,:)
    real(RP),         intent(out) :: llvely_org(:,:,:)
    real(RP),         intent(out) :: pres_org(:,:,:)
    real(RP),         intent(out) :: temp_org(:,:,:)
    real(RP),         intent(out) :: qtrc_org(:,:,:,:)
    real(RP),         intent(out) :: lon_org(:,:)
    real(RP),         intent(out) :: lat_org(:,:)
    real(RP),         intent(out) :: cz_org(:,:,:)
    character(len=*), intent(in)  :: basename
    integer,          intent(in)  :: mptype_parent
    integer,          intent(in)  :: dims(6)
    integer,          intent(in)  :: it

    ! k, i, j
    real(RP) :: velx_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: vely_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: pott_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: topo_org(          dims(2),dims(3))
    real(RP) :: geof_org(dims(4)  ,dims(2),dims(3))


    ! i, j, k
    real(RP) :: velzs_org(dims(2),dims(3),dims(4))
    real(RP) :: velxs_org(dims(5),dims(3),dims(1))
    real(RP) :: velys_org(dims(2),dims(6),dims(1))

    real(RP) :: dens
    real(RP) :: qtot

    integer :: k, i, j, iq

    character(len=H_MID) :: varname_T
    character(len=H_MID) :: varname_W
    character(len=H_MID) :: varname_U
    character(len=H_MID) :: varname_V

    integer :: fid
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputWRF]'

    if ( wrfout ) then
       varname_T = "T"
       varname_W = "W"
       varname_U = "U"
       varname_V = "V"
    else
       varname_T = "T_1"
       varname_W = "W_1"
       varname_U = "U_1"
       varname_V = "V_1"
    endif


    call FILE_open( basename, fid, rankid=myrank, single=.true., postfix="" )

    call FILE_read( fid, "XLAT", lat_org(:,:), step=it )
    lat_org(:,:) = lat_org(:,:) * D2R

    call FILE_read( fid, "XLONG", lon_org(:,:), step=it )
    lon_org(:,:) = lon_org(:,:) * D2R

    call FILE_read( fid, "HGT", topo_org(:,:), step=it )

    call FILE_read( fid, "PH", read_xyw(:,:,:), step=it )
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(4)
       ph_org(k,i,j) = read_xyw(i,j,k)
    end do
    end do
    end do

    call FILE_read( fid, "PHB", read_xyw(:,:,:), step=it )
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(4)
       phb_org(k,i,j) = read_xyw(i,j,k)
    end do
    end do
    end do

    call FILE_read( fid, "P", read_xyz(:,:,:), step=it )
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       p_org(k,i,j) = read_xyz(i,j,k)
    end do
    end do
    end do

    call FILE_read( fid, "PB", read_xyz(:,:,:), step=it )
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       pb_org(k,i,j) = read_xyz(i,j,k)
    end do
    end do
    end do

    call FILE_read( fid, varname_W, velzs_org(:,:,:), step=it )

    call FILE_read( fid, varname_U, velxs_org(:,:,:), step=it )

    call FILE_read( fid, varname_V, velys_org(:,:,:), step=it )


    ! from half level to full level
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 1, dims(1)
          velz_org(k+2,i,j) = ( velzs_org(i,j,k) + velzs_org(i,j,k+1) ) * 0.5_RP
          velx_org(k+2,i,j) = ( velxs_org(i,j,k) + velxs_org(i+1,j,k) ) * 0.5_RP
          vely_org(k+2,i,j) = ( velys_org(i,j,k) + velys_org(i,j+1,k) ) * 0.5_RP
       end do
       velz_org(1:2,i,j) = 0.0_RP
       velx_org(1:2,i,j) = 0.0_RP
       vely_org(1:2,i,j) = 0.0_RP
    end do
    end do

    call wrf_arwpost_calc_uvmet( llvelx_org, llvely_org, & ! (out)
                                 velx_org, vely_org,     & ! (in)
                                 LON_org, LAT_org,       & ! (in)
                                 BASENAME,               & ! (in)
                                 dims(1)+2, dims(2), dims(3) ) ! (in)

    qtrc_org(:,:,:,:) = 0.0_RP

    call FILE_read( fid, "Q2", read_xy(:,:), step=it )
    do j = 1, dims(3)
    do i = 1, dims(2)
       qtrc_org(1,i,j,I_QV) = read_xy(i,j)
       qtrc_org(2,i,j,I_QV) = read_xy(i,j)
    end do
    end do

    call FILE_read( fid, "QVAPOR", read_xyz(:,:,:), step=it )
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       qtrc_org(k+2,i,j,I_QV) = read_xyz(i,j,k)
    end do
    end do
    end do


    if ( I_QC > 0 ) then
       call FILE_read( fid, "QCLOUD", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_QC) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    if ( I_QR > 0 ) then
       call FILE_read( fid, "QRAIN", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_QC) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    if ( I_QI > 0 ) then
       call FILE_read( fid, "QICE", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_QI) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    if ( I_QS > 0 ) then
       call FILE_read( fid, "QSNOW", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_QS) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    if ( I_QG > 0 ) then
       call FILE_read( fid, "QGRAUP", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_QG) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if


    ! convert mixing ratio to specific ratio
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)+2
       qtot = qtrc_org(k,i,j,I_QV)
       do iq = QHS, QHE
          qtot = qtot + qtrc_org(k,i,j,iq)
       end do
       qtrc_org(k,i,j,I_QV) = qtrc_org(k,i,j,I_QV) / ( 1.0_RP + qtot )
       do iq = QHS, QHE
          qtrc_org(k,i,j,iq) = qtrc_org(k,i,j,iq) / ( 1.0_RP + qtot )
       end do
    end do
    end do
    end do

    if ( I_NC > 0 ) then
       call FILE_read( fid, "NC", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_NC) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    if ( I_NR > 0 ) then
       call FILE_read( fid, "NR", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_NR) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    if ( I_NI > 0 ) then
       call FILE_read( fid, "NI", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_NI) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    if ( I_NS > 0 ) then
       call FILE_read( fid, "NS", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_NS) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    if ( I_NG > 0 ) then
       call FILE_read( fid, "NG", read_xyz(:,:,:), step=it, allow_missing=.true. )
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          qtrc_org(k+2,i,j,I_NG) = read_xyz(i,j,k)
       end do
       end do
       end do
    end if

    do iq = 1, QA
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)+2
       qtrc_org(k,i,j,iq) = max( qtrc_org(k,i,j,iq), 0.0_RP )
    end do
    end do
    end do
    end do


    call FILE_read( fid, varname_T, read_xyz(:,:,:), step=it, allow_missing=.true. )
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       pott_org(k+2,i,j) = read_xyz(i,j,k) + t0
    end do
    end do
    end do
    call FILE_read( fid, "T2", read_xy(:,:), step=it, allow_missing=.true. )
    do j = 1, dims(3)
    do i = 1, dims(2)
       temp_org(2,i,j) = read_xy(i,j)
    end do
    end do

    call FILE_read( fid, "PSFC", read_xy(:,:), step=it, allow_missing=.true. )
    do j = 1, dims(3)
    do i = 1, dims(2)
       pres_org(2,i,j) = read_xy(i,j)
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 3, dims(1)+2
          pres_org(k,i,j) = p_org(k-2,i,j) + pb_org(k-2,i,j)
          temp_org(k,i,j) = pott_org(k,i,j) * ( pres_org(k,i,j) / p0 )**RCP
       end do
       pott_org(2,i,j) = temp_org(2,i,j) * ( p0/pres_org(2,i,j) )**RCP
       temp_org(1,i,j) = temp_org(2,i,j) + LAPS * topo_org(i,j)
       dens = pres_org(2,i,j) / ( Rdry * temp_org(2,i,j) )
       pres_org(1,i,j) = ( pres_org(2,i,j) + GRAV * dens * cz_org(2,i,j) * 0.5_RP ) &
                       / ( Rdry * temp_org(1,i,j) - GRAV * cz_org(2,i,j) * 0.5_RP ) &
                       * Rdry * temp_org(1,i,j)
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)
       ! convert to geopotential height to use as real height in WRF
       do k = 1, dims(4)
          geof_org(k,i,j) = ( ph_org(k,i,j) + phb_org(k,i,j) ) / grav
       end do
       ! make half level of geopotential height from face level
       do k = 1, dims(1)
          cz_org(k+2,i,j) = ( geof_org(k,i,j) + geof_org(k+1,i,j) ) * 0.5_RP
       end do
       cz_org(2,i,j) = topo_org(i,j)
       cz_org(1,i,j) = 0.0_RP
    end do
    end do


#ifdef DEBUG
     !k=1 ; i=int(dims(2)/2) ; j=int(dims(3)/2) ; iq = 1
     k=2 ; i=3 ; j=3 ; iq = 1
     write(*,*) "read 3D wrf data",i,j,k
     write(*,*) "lon_org    ",lon_org   (i,j)/D2R
     write(*,*) "lat_org    ",lat_org   (i,j)/D2R
     write(*,*) "cz_org     ",cz_org    (k,i,j)
     write(*,*) "pres_org   ",pres_org  (k,i,j)
     write(*,*) "velx_org   ",llvelx_org(k,i,j)
     write(*,*) "vely_org   ",llvely_org(k,i,j)
     write(*,*) "velz_org   ",velz_org  (k,i,j)
     write(*,*) "temp_org   ",temp_org  (k,i,j)
     write(*,*) "qtrc_org   ",qtrc_org  (k,i,j,iq)
     k=3 ; i=3 ; j=3 ; iq = 1
     write(*,*) "read 3D wrf data",i,j,k
     write(*,*) "lon_org    ",lon_org   (i,j)/D2R
     write(*,*) "lat_org    ",lat_org   (i,j)/D2R
     write(*,*) "cz_org     ",cz_org    (k,i,j)
     write(*,*) "pres_org   ",pres_org  (k,i,j)
     write(*,*) "velx_org   ",llvelx_org(k,i,j)
     write(*,*) "vely_org   ",llvely_org(k,i,j)
     write(*,*) "velz_org   ",velz_org  (k,i,j)
     write(*,*) "temp_org   ",temp_org  (k,i,j)
     write(*,*) "qtrc_org   ",qtrc_org  (k,i,j,iq)
#endif

    return
  end subroutine ParentAtmosInputWRFARW

  !-----------------------------------------------------------------------------
  !> Land Setup
  subroutine ParentLandSetupWRFARW( &
       ldims,    &
       basename_land )
    use scale_file, only: &
       FILE_open, &
       FILE_get_dimLength
    implicit none

    integer,          intent(out) :: ldims(3)
    character(len=*), intent(in)  :: basename_land

    logical :: WRF_FILE_TYPE = .false.   ! wrf filetype: T=wrfout, F=wrfrst

    NAMELIST / PARAM_MKINIT_REAL_WRFARW / &
         WRF_FILE_TYPE

    integer :: fid
    integer :: timelen
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atmos Input File Type: WRF-ARW'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_WRFARW,iostat=ierr)
    if( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_WRFARW. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_MKINIT_REAL_WRFARW)


    call FILE_open( basename_land, fid, rankid=myrank, single=.true., postfix="" )

    call FILE_get_dimLength( fid, "soil_layers_stag", ldims(1) )
    call FILE_get_dimLength( fid, "west_east",        ldims(2) )
    call FILE_get_dimLength( fid, "south_north",      ldims(3) )

    if ( wrf_file_type ) then
       wrfout = .true.
       if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF History Output'
    else
       wrfout = .false.
       if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF Restart'
    endif


    if ( .not. allocated(read_xy) ) then
       allocate( read_xy  (ldims(2),ldims(3)) )
    end if

    allocate( read_xyl(ldims(2),ldims(3),ldims(1)) )

    return
  end subroutine ParentLandSetupWRFARW

  !-----------------------------------------------------------------------------
  subroutine ParentLandInputWRFARW( &
      tg_org,             &
      sh2o_org,           &
      lst_org,            &
      ust_org,            &
      albg_org,           &
      topo_org,           &
      lmask_org,          &
      llon_org,           &
      llat_org,           &
      lz_org,             &
      basename,           &
      ldims,              &
      use_file_landwater, &
      it                  )
    use scale_const, only: &
       D2R => CONST_D2R, &
       UNDEF => CONST_UNDEF, &
       I_LW => CONST_I_LW, &
       I_SW => CONST_I_SW
    use scale_file, only: &
       FILE_open, &
       FILE_read
    implicit none
    real(RP),         intent(out)  :: tg_org(:,:,:)
    real(RP),         intent(out)  :: sh2o_org(:,:,:)
    real(RP),         intent(out)  :: lst_org(:,:)
    real(RP),         intent(out)  :: ust_org(:,:)
    real(RP),         intent(out)  :: albg_org(:,:,:)
    real(RP),         intent(out)  :: topo_org(:,:)
    real(RP),         intent(out)  :: lmask_org(:,:)
    real(RP),         intent(out)  :: llon_org(:,:)
    real(RP),         intent(out)  :: llat_org(:,:)
    real(RP),         intent(out)  :: lz_org(:)
    character(len=*), intent( in)  :: basename
    integer,          intent( in)  :: ldims(3)
    logical,          intent( in)  :: use_file_landwater   ! use land water data from files
    integer,          intent( in)  :: it

    integer :: fid
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[LandInputWRF]'

    call FILE_open( basename, fid, rankid=myrank, single=.true., postfix="" )

    call FILE_read( fid, "XLAT", llat_org(:,:), step=it )
    llat_org(:,:) = llat_org(:,:) * D2R

    call FILE_read( fid, "XLONG", llon_org(:,:), step=it )
    llon_org(:,:) = llon_org(:,:) * D2R

    call FILE_read( fid, "HGT", topo_org(:,:), step=it )


    ! depth
    call FILE_read( fid, "ZS", lz_org(:), step=it )

    ! land mask (1:land, 0:water)
    call FILE_read( fid, "LANDMASK", lmask_org(:,:), step=it )

    ! soil temperature [K]
    call FILE_read( fid, "TSLB", read_xyl(:,:,:), step=it )
    do j = 1, ldims(3)
    do i = 1, ldims(2)
    do k = 1, ldims(1)
       tg_org(k,i,j) = read_xyl(i,j,k)
    end do
    end do
    end do

    ! soil liquid water [m3 m-3] (no wrfout-default)
    if( use_file_landwater ) then
       call FILE_read( fid, "SH2O", read_xyl(:,:,:), step=it, allow_missing=.true., missing_value=UNDEF )
       do j = 1, ldims(3)
       do i = 1, ldims(2)
       do k = 1, ldims(1)
          sh2o_org(k,i,j) = read_xyl(i,j,k)
       end do
       end do
       end do
    endif

!    ! surface runoff [mm]
!    call FILE_read( fid, "SFROFF", org_3D(:,:), step=it )
!    do j = 1, ldims(3)
!    do i = 1, ldims(2)
!       org_3D(k,i,j) = org_3D(i,j,k) * 1000.0_RP * dwatr
!    end do
!    end do


    ! SURFACE SKIN TEMPERATURE [K]
    call FILE_read( fid, "TSK", lst_org(:,:), step=it )

    ust_org(:,:) = lst_org(:,:)

    ! ALBEDO [-]
    call FILE_read( fid, "ALBEDO", albg_org(:,:,I_SW), step=it )

    ! SURFACE EMISSIVITY [-]
    call FILE_read( fid, "EMISS", read_xy(:,:), step=it )
    do j = 1, ldims(3)
    do i = 1, ldims(2)
       albg_org(i,j,I_LW) = 1.0_RP - read_xy(i,j)
    end do
    end do


!    ! SNOW WATER EQUIVALENT [kg m-2] (no wrfout-default)
!    call FILE_read( fid, "SNOW", snowq_org(:,:), step=it, allow_missing=.true., missing_value=UNDEF )

!    ! AVERAGE SNOW TEMPERATURE [C] (no wrfout-default)
!    call FILE_read( fid, "TSNAV", snowt_org(:,:), step=it, allow_missing=.true., missing_value=UNDEF )
!    do j = 1, ldims(3)
!    do i = 1, ldims(2)
!       if ( snowt_org(k,i,j) /= UNDEF ) snowt_org(k,i,j) = snowt_org(i,j,k) + TEM00
!    end do
!    end do

    return
  end subroutine ParentLandInputWRFARW

  !-----------------------------------------------------------------------------
  !> Ocean Setup
  subroutine ParentOceanSetupWRFARW( &
       odims,    &
       timelen, &
       basename_org )
    use scale_file, only: &
       FILE_open, &
       FILE_get_dimLength
    implicit none

    integer,          intent(out) :: odims(2)
    integer,          intent(out) :: timelen
    character(len=*), intent(in)  :: basename_org

    logical :: WRF_FILE_TYPE = .false.   ! wrf filetype: T=wrfout, F=wrfrst

    NAMELIST / PARAM_MKINIT_REAL_WRFARW / &
         WRF_FILE_TYPE

    integer :: fid
    integer :: ierr
    logical :: error
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Ocean Input File Type: WRF-ARW'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_WRFARW,iostat=ierr)
    if( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_WRFARW. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_MKINIT_REAL_WRFARW)


    call FILE_open( basename_org, fid, rankid=myrank, single=.true., postfix="" )

    call FILE_get_dimLength(fid, "west_east",   odims(1) )
    call FILE_get_dimLength(fid, "south_north", odims(2) )

    call FILE_get_dimLength(fid, "Time", timelen, error=error )
    if ( error ) call FILE_get_dimLength(fid, "time",  timelen, error=error)
    if ( error ) timelen = 0

    if ( wrf_file_type ) then
       wrfout = .true.
       if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF History Output'
    else
       wrfout = .false.
       if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF Restart'
    endif


    if ( .not. allocated(read_xy) ) then
       allocate( read_xy(odims(1),odims(2)) )
    end if

    return
  end subroutine ParentOceanSetupWRFARW

  !-----------------------------------------------------------------------------
  subroutine ParentOceanOpenWRFARW
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OceanOpenWRFARW]'
    return
  end subroutine ParentOceanOpenWRFARW

  !-----------------------------------------------------------------------------
  subroutine ParentOceanInputWRFARW( &
      tw_org,             &
      sst_org,            &
      albw_org,           &
      z0w_org,            &
      omask_org,          &
      olon_org,           &
      olat_org,           &
      basename,           &
      odims,              &
      it                  )
    use scale_const, only: &
       D2R => CONST_D2R, &
       UNDEF => CONST_UNDEF, &
       I_LW => CONST_I_LW, &
       I_SW => CONST_I_SW
    use scale_file, only: &
       FILE_open, &
       FILE_read
    implicit none
    real(RP),         intent(out)  :: tw_org(:,:)
    real(RP),         intent(out)  :: sst_org(:,:)
    real(RP),         intent(out)  :: albw_org(:,:,:)
    real(RP),         intent(out)  :: z0w_org(:,:)
    real(RP),         intent(out)  :: omask_org(:,:)
    real(RP),         intent(out)  :: olon_org(:,:)
    real(RP),         intent(out)  :: olat_org(:,:)
    character(len=*), intent( in)  :: basename
    integer,          intent( in)  :: odims(2)
    integer,          intent( in)  :: it

    integer :: fid
    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OceanInputWRF]'

    call FILE_open( basename, fid, rankid=myrank, single=.true., postfix="" )

    call FILE_read( fid, "XLAT", olat_org(:,:), step=it )
    olat_org(:,:) = olat_org(:,:) * D2R

    call FILE_read( fid, "XLONG", olon_org(:,:), step=it )
    olon_org(:,:) = olon_org(:,:) * D2R


    ! land mask (1:land, 0:water)
    call FILE_read( fid, "LANDMASK", omask_org(:,:), step=it )

    ! SEA SURFACE TEMPERATURE [K]
    call FILE_read( fid, "SST", sst_org(:,:), step=it )

    tw_org(:,:) = sst_org(:,:)

    ! ALBEDO [-]
    call FILE_read( fid, "ALBEDO", albw_org(:,:,I_SW), step=it )

    ! SURFACE EMISSIVITY [-]
    call FILE_read( fid, "EMISS", read_xy(:,:), step=it )
    do j = 1, odims(2)
    do i = 1, odims(1)
       albw_org(i,j,I_LW) = 1.0_RP - read_xy(i,j)
    enddo
    enddo

    ! TIME-VARYING ROUGHNESS LENGTH [m] (no wrfout-default)
    call FILE_read( fid, "ZNT", z0w_org(:,:), step=it, allow_missing=.true., missing_value=UNDEF )


    return
  end subroutine ParentOceanInputWRFARW

  !-----------------------------------------------------------------------------
  !> convert vector varibles from map-projected grid on wrf model to lat-lon grid
  !-----------------------------------------------------------------------------
  subroutine wrf_arwpost_calc_uvmet( &
      u_latlon,   & ! (out)
      v_latlon,   & ! (out)
      u_on_map,   & ! (in)
      v_on_map,   & ! (in)
      xlon,       & ! (in)
      xlat,       & ! (in)
      basename ,  & ! (in)
      K1, I1, J1  ) ! (in)
    use scale_const, only: &
       D2R => CONST_D2R, &
       PI => CONST_PI
    use scale_file, only: &
       FILE_open, &
       FILE_get_attribute
    implicit none
    real(RP), intent(out) :: u_latlon(:,:,:)
    real(RP), intent(out) :: v_latlon(:,:,:)
    real(RP), intent(in ) :: u_on_map(:,:,:)
    real(RP), intent(in ) :: v_on_map(:,:,:)
    real(RP), intent(in ) :: xlon(:,:)
    real(RP), intent(in ) :: xlat(:,:)
    integer , intent(in ) :: K1, I1, J1

    character(len=*), intent( in) :: basename

    integer :: fid

    real(RP) :: truelat1, truelat2
    real(RP) :: stand_lon
    real(RP) :: diff
    real(RP) :: alpha
    real(RP) :: sine(I1,J1)
    real(RP) :: cose(I1,J1)
    real(RP) :: cone
    integer  :: map_proj

    real(RP) :: dum_r(1)
    integer  :: dum_i(1)


    integer  :: k, i, j
    !---------------------------------------------------------------------------

    call FILE_open( basename, fid, rankid=myrank, single=.true., postfix="" )

    call FILE_get_attribute( fid, "global", "MAP_PROJ", dum_i(:) )
    map_proj = dum_i(1)

    call FILE_get_attribute( fid, "global", "TRUELAT1", dum_r(:) )
    truelat1 = dum_r(1) * D2R
    call FILE_get_attribute( fid, "global", "TRUELAT2", dum_r(:) )
    truelat2 = dum_r(1) * D2R
    call FILE_get_attribute( fid, "global", "STAND_LON", dum_r(:) )
    stand_lon = dum_r(1) * D2R

    ! No need to rotate
    if ( map_proj .ge. 3 ) then
       u_latlon(:,:,:) = u_on_map(:,:,:)
       v_latlon(:,:,:) = v_on_map(:,:,:)

       return
    endif

    ! Lambert Conformal mapping
    cone = 1.0_RP         !  PS
    if ( map_proj .eq. 1 ) then
       if ( abs(truelat1-truelat2) .gt. 0.1_RP*D2R ) then
          cone = ( log(cos(truelat1)) -                 &
                   log(cos(truelat2)) )  /              &
                 ( log(tan((PI*0.5_RP-abs(truelat1))*0.5_RP )) - &
                   log(tan((PI*0.5_RP-abs(truelat2))*0.5_RP )) )
       else
          cone = sin( abs(truelat1) )
       endif
    endif

    do j = 1, J1
    do i = 1, I1
       diff = xlon(i,j) - stand_lon
       if ( diff .gt. PI ) then
          diff = diff - PI*2.0_RP
       endif
       if ( diff .lt. -PI ) then
          diff = diff + PI*2.0_RP
       endif
       alpha = diff * cone * sign(1.0_RP, xlat(i,j))
       sine(i,j) = sin( alpha )
       cose(i,j) = cos( alpha )
    enddo
    enddo

    do j = 1, J1
    do i = 1, I1
    do k = 1, K1
       u_latlon(k,i,j) = v_on_map(k,i,j)*sine(i,j) + u_on_map(k,i,j)*cose(i,j)
       v_latlon(k,i,j) = v_on_map(k,i,j)*cose(i,j) - u_on_map(k,i,j)*sine(i,j)
    enddo
    enddo
    enddo

    return
  end subroutine wrf_arwpost_calc_uvmet

end module mod_realinput_wrfarw
