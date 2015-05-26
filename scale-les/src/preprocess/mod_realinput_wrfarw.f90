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
     PRC_master,            &
     PRC_MPIstop
  use scale_external_io, only: &
     iWRFARW

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtomSetupWRFARW
  public :: ParentAtomOpenWRFARW
  public :: ParentAtomInputWRFARW
  public :: ParentSurfaceInputWRFARW

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
  integer,  parameter :: mdlid = iWRFARW

  real(RP), parameter :: t0  = 300.0_RP
  real(RP), parameter :: p0  = 1000.0E+2_RP
  real(RP), parameter :: Rd  = 287.04_RP
  real(RP), parameter :: Cp  = 7.0_RP * Rd / 2.0_RP
  real(RP), parameter :: RCP = Rd / Cp

  integer, parameter :: cosin = 1
  integer, parameter :: sine  = 2

  real(SP), allocatable :: read_xy (:,:,:)
  real(SP), allocatable :: read_uy (:,:,:)
  real(SP), allocatable :: read_xv (:,:,:)
  real(SP), allocatable :: read_zxy(:,:,:,:)
  real(SP), allocatable :: read_wxy(:,:,:,:)
  real(SP), allocatable :: read_zuy(:,:,:,:)
  real(SP), allocatable :: read_zxv(:,:,:,:)
  real(SP), allocatable :: read_lzxy(:,:,:,:)
  real(SP), allocatable :: read_lz(:,:)

  real(RP), allocatable :: p_org   (:,:,:)
  real(RP), allocatable :: pb_org  (:,:,:)
  real(RP), allocatable :: ph_org  (:,:,:)
  real(RP), allocatable :: phb_org (:,:,:)

  logical, private :: wrfout = .false. ! file type switch (wrfout or wrfrst)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ParentAtomSetupWRFARW( &
      dims,    &
      timelen, &
      basename_org )
    use scale_external_io, only: &
         iWRFARW, &
         ExternalFileGetShape
    implicit none

    integer,               intent(out) :: dims(11)
    integer,               intent(out) :: timelen
    character(len=H_LONG), intent(in) :: basename_org

    logical :: WRF_FILE_TYPE = .false.   ! wrf filetype: T=wrfout, F=wrfrst

    NAMELIST / PARAM_MKINIT_REAL_WRFARW / &
         WRF_FILE_TYPE

    integer :: dims_wrf(7)
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: WRF-ARW'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_WRFARW,iostat=ierr)
    if( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_WRFARW. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_REAL_WRFARW)


    call ExternalFileGetShape( dims_wrf, timelen, mdlid, basename_org, myrank, single=.true. )
    dims(1:7) = dims_wrf
    ! land
    dims(8)  = dims(2)
    dims(9)  = dims(3)
    ! sst
    dims(10) = dims(2)
    dims(11) = dims(3)

    if ( wrf_file_type ) then
       wrfout = .true.
       if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF History Output'
    else
       wrfout = .false.
       if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF Restart'
    endif


    allocate( read_xy  (        dims(2),dims(3),1) )
    allocate( read_uy  (        dims(5),dims(3),1) )
    allocate( read_xv  (        dims(2),dims(6),1) )
    allocate( read_zxy (dims(1),dims(2),dims(3),1) )
    allocate( read_wxy (dims(4),dims(2),dims(3),1) )
    allocate( read_zuy (dims(1),dims(5),dims(3),1) )
    allocate( read_zxv (dims(1),dims(2),dims(6),1) )

    allocate( read_lzxy(dims(7),dims(8),dims(9),1) )
    allocate( read_lz  (dims(7),1) )

    allocate( p_org    (dims(1),dims(2),dims(3)) )
    allocate( pb_org   (dims(1),dims(2),dims(3)) )
    allocate( ph_org   (dims(4),dims(2),dims(3)) )
    allocate( phb_org  (dims(4),dims(2),dims(3)) )

    return
  end subroutine ParentAtomSetupWRFARW

  !-----------------------------------------------------------------------------
  subroutine ParentAtomOpenWRFARW
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OpenWRFARW]'
    return
  end subroutine ParentAtomOpenWRFARW

  !-----------------------------------------------------------------------------
  subroutine ParentAtomInputWRFARW( &
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
         GRAV => CONST_GRAV
    use scale_external_io, only: &
         ExternalFileRead
    use scale_atmos_thermodyn, only: &
         THERMODYN_pott => ATMOS_THERMODYN_pott
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
    character(LEN=*), intent(in)  :: basename
    integer,          intent(in)  :: mptype_parent
    integer,          intent(in)  :: dims(7)
    integer,          intent(in)  :: it

    ! full level
    real(RP) :: velx_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: vely_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: topo_org(          dims(2),dims(3))

    ! half level
    real(RP) :: velzs_org(dims(4),dims(2),dims(3))
    real(RP) :: velxs_org(dims(1),dims(5),dims(3))
    real(RP) :: velys_org(dims(1),dims(2),dims(6))
    real(RP) :: geof_org (dims(4),dims(2),dims(3))

    real(RP) :: pott
    real(RP) :: qhyd

    integer :: k, i, j, iq

    character(len=H_MID) :: varname_T
    character(len=H_MID) :: varname_W
    character(len=H_MID) :: varname_U
    character(len=H_MID) :: varname_V

    logical :: lack_of_val

    integer :: ierr
    !---------------------------------------------------------------------------

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


    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputWRF]'

    call ExternalFileRead( read_xy(:,:,:),    BASENAME, "XLAT",    it, it, myrank, mdlid, single=.true.               )
    lat_org (:,:) = real( read_xy(:,:,1), kind=RP ) * D2R

    call ExternalFileRead( read_xy(:,:,:),    BASENAME, "XLONG",   it, it, myrank, mdlid, single=.true.               )
    lon_org (:,:) = real( read_xy(:,:,1), kind=RP ) * D2R

    call ExternalFileRead( read_xy(:,:,:),    BASENAME, "HGT",     it, it, myrank, mdlid, single=.true.               )
    topo_org (:,:) = real( read_xy(:,:,1), kind=RP )

    call ExternalFileRead( read_wxy(:,:,:,:), BASENAME, "PH",      it, it, myrank, mdlid, single=.true., zstag=.true. )
    ph_org  (3:,:,:) = real( read_wxy(:,:,:,1), kind=RP )

    call ExternalFileRead( read_wxy(:,:,:,:), BASENAME, "PHB",     it, it, myrank, mdlid, single=.true., zstag=.true. )
    phb_org (3:,:,:) = real( read_wxy(:,:,:,1), kind=RP )

    call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "P",       it, it, myrank, mdlid, single=.true.               )
    p_org   (3:,:,:) = real( read_zxy(:,:,:,1), kind=RP )

    call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "PB",      it, it, myrank, mdlid, single=.true.               )
    pb_org  (3:,:,:) = real( read_zxy(:,:,:,1), kind=RP )

    call ExternalFileRead( read_wxy(:,:,:,:), BASENAME, varname_W, it, it, myrank, mdlid, single=.true., zstag=.true. )
    velzs_org(:,:,:) = real( read_wxy(:,:,:,1), kind=RP )

    call ExternalFileRead( read_zuy(:,:,:,:), BASENAME, varname_U, it, it, myrank, mdlid, single=.true., xstag=.true. )

    velxs_org(:,:,:) = real( read_zuy(:,:,:,1), kind=RP )

    call ExternalFileRead( read_zxv(:,:,:,:), BASENAME, varname_V, it, it, myrank, mdlid, single=.true., ystag=.true. )
    velys_org(:,:,:) = real( read_zxv(:,:,:,1), kind=RP )

    ! from half level to full level
    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 1, dims(1)
          velz_org(k+2,i,j) = ( velzs_org(k,i,j) + velzs_org(k+1,i,j) ) * 0.5_RP
          velx_org(k+2,i,j) = ( velxs_org(k,i,j) + velxs_org(k,i+1,j) ) * 0.5_RP
          vely_org(k+2,i,j) = ( velys_org(k,i,j) + velys_org(k,i,j+1) ) * 0.5_RP
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
    call ExternalFileRead( read_xy(:,:,:),    BASENAME, "Q2",      it, it, myrank, mdlid, single=.true.               )
    qtrc_org(2,:,:,I_QV) = real(read_xy(:,:,1),kind=RP)
    qtrc_org(1,:,:,I_QV) = qtrc_org(2,:,:,I_QV)

    call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QVAPOR",  it, it, myrank, mdlid, single=.true. )
    qtrc_org(3:,:,:,I_QV) = real(read_zxy(:,:,:,1),kind=RP)

#ifndef DRY
    if( mptype_parent > 0 ) then
       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QCLOUD", it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_QC) = real( read_zxy(:,:,:,1), kind=RP)

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QRAIN",  it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_QR) = real( read_zxy(:,:,:,1), kind=RP)
    endif

    if( mptype_parent > 3 ) then
       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QICE",   it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_QI) = real( read_zxy(:,:,:,1), kind=RP)

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QSNOW",  it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_QS) = real( read_zxy(:,:,:,1), kind=RP)
    endif

    if( mptype_parent > 5 ) then
       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QGRAUP", it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_QG) = real( read_zxy(:,:,:,1), kind=RP)
    endif

    ! convert mixing ratio to specific ratio
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)+2
       qhyd = 0.0_RP
       do iq = 1, min( mptype_parent, 6)
          qhyd = qhyd + qtrc_org(k,i,j,iq)
       end do
       do iq = 1, min( mptype_parent, 6)
          qtrc_org(k,i,j,iq) = qtrc_org(k,i,j,iq) * ( 1.0_RP - qhyd )
       end do
    end do
    end do
    end do

    if( mptype_parent > 6 ) then
       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NC",     it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_NC) = real( read_zxy(:,:,:,1), kind=RP )

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NR",     it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_NR) = real( read_zxy(:,:,:,1), kind=RP )

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NI",     it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_NI) = real( read_zxy(:,:,:,1), kind=RP )

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NS",     it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_NS) = real( read_zxy(:,:,:,1), kind=RP )

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NG",     it, it, myrank, mdlid, single=.true. )
       qtrc_org(3:,:,:,I_NG) = real( read_zxy(:,:,:,1), kind=RP )
    endif
#endif

    do iq = 1, QA
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)+2
       qtrc_org(k,i,j,iq) = max( qtrc_org(k,i,j,iq), 0.0_RP )
    end do
    end do
    end do
    end do


    call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, varname_T, it, it, myrank, mdlid, single=.true.               )
    temp_org(3:,:,:) = real( read_zxy(:,:,:,1), kind=RP )
    call ExternalFileRead( read_xy(:,:,:),    BASENAME, "T2",      it, it, myrank, mdlid, single=.true.               )
    temp_org(2,:,:) = real( read_xy(:,:,1), kind=RP )

    call ExternalFileRead( read_xy(:,:,:),    BASENAME, "PSFC",    it, it, myrank, mdlid, single=.true.               )
    pres_org(2,:,:) = real( read_xy(:,:,1), kind=RP )

    do j = 1, dims(3)
    do i = 1, dims(2)
       do k = 3, dims(1)+2
          pres_org(k,i,j) = p_org(k,i,j) + pb_org(k,i,j)
       end do
       pott = temp_org(2,i,j) * ( p0/pres_org(2,i,j) )**RCP
       temp_org(1,i,j) = temp_org(2,i,j) + LAPS * topo_org(i,j)
       pres_org(1,i,j) = p0 * ( temp_org(1,i,j) / pott )**(1.0_RP/RCP)
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


    return
  end subroutine ParentAtomInputWRFARW

  !-----------------------------------------------------------------------------
  subroutine ParentSurfaceInputWRFARW( &
      tg_org,             &
      sh2o_org,           &
      tw_org,             &
      lst_org,            &
      ust_org,            &
      sst_org,            &
      albw_org,           &
      albg_org,           &
      z0w_org,            &
      lmask_org,          &
      lz_org,             &
      basename,           &
      dims,               &
      use_file_landwater, &
      it                  )
    use scale_const, only: &
         UNDEF => CONST_UNDEF, &
         I_LW => CONST_I_LW, &
         I_SW => CONST_I_SW
    use scale_external_io, only: &
         ExternalFileRead, &
         ExternalFileVarExistence
    implicit none
    real(RP),         intent(out)  :: tg_org(:,:,:)
    real(RP),         intent(out)  :: sh2o_org(:,:,:)
    real(RP),         intent(out)  :: tw_org(:,:)
    real(RP),         intent(out)  :: lst_org(:,:)
    real(RP),         intent(out)  :: ust_org(:,:)
    real(RP),         intent(out)  :: sst_org(:,:)
    real(RP),         intent(out)  :: albw_org(:,:,:)
    real(RP),         intent(out)  :: albg_org(:,:,:)
    real(RP),         intent(out)  :: z0w_org(:,:)
    real(RP),         intent(out)  :: lmask_org(:,:)
    real(RP),         intent(out)  :: lz_org(:)
    character(LEN=*), intent( in)  :: basename
    integer,          intent( in)  :: dims(11)
    logical,          intent( in)  :: use_file_landwater   ! use land water data from files
    integer,          intent( in)  :: it


    integer  :: k, i, j, iq, iqw

    logical  :: existence

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputWRF-Surface]'

    ! depth
    call ExternalFileRead( read_lz(:,:),                               &
                      BASENAME, "ZS",      it, 1, myrank, mdlid, dims(7), single=.true. )
    lz_org(:) = read_lz(:,1)

    ! land mask (1:land, 0:water)
    call ExternalFileRead( read_xy(:,:,:),                             &
                      BASENAME, "LANDMASK",  it, 1, myrank, mdlid, single=.true. )
    lmask_org(:,:) = read_xy(:,:,1)

    ! soil temperature [K]
    call ExternalFileRead( read_lzxy(:,:,:,:),                             &
                      BASENAME, "TSLB",  it, 1, myrank, mdlid, single=.true., landgrid=.true. )
    tg_org(:,:,:) = read_lzxy(:,:,:,1)

    ! soil liquid water [m3 m-3] (no wrfout-default)
    if( use_file_landwater ) then
       call ExternalFileVarExistence( existence, BASENAME, "SH2O", myrank, mdlid, single=.true. )
       if ( existence ) then
          call ExternalFileRead( read_lzxy(:,:,:,:),                             &
                         BASENAME, "SH2O", it, 1, myrank, mdlid, single=.true., landgrid=.true.  )
          sh2o_org(:,:,:) = read_lzxy(:,:,:,1)
       else
          sh2o_org(:,:,:) = UNDEF
       endif
    endif

!    ! surface runoff [mm]
!    call ExternalFileRead( read_xy(:,:,:),                             &
!                      BASENAME, "SFROFF",  it, 1, myrank, mdlid, single=.true. )
!       do j = 1, dims(9)
!       do i = 1, dims(8)
!          org_3D(i,j) = read_xy(i,j,1) * 1000.0_DP * dwatr
!       enddo
!       enddo


    ! SEA SURFACE TEMPERATURE [K]
    call ExternalFileRead( read_xy(:,:,:),                             &
                      BASENAME, "SST",  it, 1, myrank, mdlid, single=.true. )
    tw_org(:,:) = read_xy(:,:,1)

    sst_org(:,:) = tw_org(:,:)

    ! SURFACE SKIN TEMPERATURE [K]
    call ExternalFileRead( read_xy(:,:,:),                             &
                      BASENAME, "TSK",  it, 1, myrank, mdlid, single=.true. )
    lst_org(:,:) = read_xy(:,:,1)

    ust_org(:,:) = lst_org(:,:)

    ! ALBEDO [-]
    call ExternalFileRead( read_xy(:,:,:),                             &
                      BASENAME, "ALBEDO",  it, 1, myrank, mdlid, single=.true. )
    albw_org(:,:,I_SW) = read_xy(:,:,1)
    albg_org(:,:,I_SW) = albw_org(:,:,I_SW)

    ! SURFACE EMISSIVITY [-]
    call ExternalFileRead( read_xy(:,:,:),                             &
                      BASENAME, "EMISS",  it, 1, myrank, mdlid, single=.true. )
    do j = 1, dims(9)
    do i = 1, dims(8)
       albw_org(i,j,I_LW) = 1.0_DP - read_xy(i,j,1)
       albg_org(i,j,I_LW) = 1.0_DP - read_xy(i,j,1)
    enddo
    enddo

    ! TIME-VARYING ROUGHNESS LENGTH [m] (no wrfout-default)
    call ExternalFileVarExistence( existence, BASENAME, "ZNT", myrank, mdlid, single=.true. )
    if ( existence ) then
       call ExternalFileRead( read_xy(:,:,:),                             &
                         BASENAME, "ZNT",  it, 1, myrank, mdlid, single=.true. )
       z0w_org(:,:) = read_xy(:,:,1)
    else
       z0w_org(:,:) = UNDEF
    endif


!    ! SNOW WATER EQUIVALENT [kg m-2] (no wrfout-default)
!    call ExternalFileVarExistence( existence, BASENAME, "SNOW", myrank, mdlid, single=.true. )
!    if ( existence ) then
!       call ExternalFileRead( read_xy(:,:,:),                             &
!                         BASENAME, "SNOW",  it, 1, myrank, mdlid, single=.true. )
!       swowq_org(:,:,:) = read_xy(:,:,1)
!    else
!       snowq_org(:,:) = UNDEF
!    endif

!    ! AVERAGE SNOW TEMPERATURE [C] (no wrfout-default)
!    call ExternalFileVarExistence( existence, BASENAME, "TSNAV", myrank, mdlid, single=.true. )
!    if ( existence ) then
!       call ExternalFileRead( read_xy(:,:,:),                             &
!                         BASENAME, "TSNAV",  it, 1, myrank, mdlid, single=.true. )
!       do j = 1, dims(9)
!       do i = 1, dims(8)
!          snowt_org(i,j,n) = read_xy(i,j,1) + TEM00
!       enddo
!       enddo
!    else
!       snowt_org(:,:) = UNDEFF
!    endif


    return
  end subroutine ParentSurfaceInputWRFARW

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
    use scale_external_io, only: &
         ExternalFileGetGlobalAttV
    implicit none
    real(RP), intent(out) :: u_latlon(:,:,:)
    real(RP), intent(out) :: v_latlon(:,:,:)
    real(RP), intent(in ) :: u_on_map(:,:,:)
    real(RP), intent(in ) :: v_on_map(:,:,:)
    real(RP), intent(in ) :: xlon(:,:)
    real(RP), intent(in ) :: xlat(:,:)
    integer , intent(in ) :: K1, I1, J1

    character(LEN=*), intent( in) :: basename

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

    call ExternalFileGetGlobalAttV( dum_i, iWRFARW, BASENAME, "MAP_PROJ",  myrank, single=.true. )
    map_proj = dum_i(1)
    call ExternalFileGetGlobalAttV( dum_r, iWRFARW, BASENAME, "TRUELAT1",  myrank, single=.true. )
    truelat1 = dum_r(1) * D2R
    call ExternalFileGetGlobalAttV( dum_r, iWRFARW, BASENAME, "TRUELAT2",  myrank, single=.true. )
    truelat2 = dum_r(1) * D2R
    call ExternalFileGetGlobalAttV( dum_r, iWRFARW, BASENAME, "STAND_LON", myrank, single=.true. )
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
!-------------------------------------------------------------------------------
