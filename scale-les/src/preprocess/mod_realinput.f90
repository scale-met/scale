!-------------------------------------------------------------------------------
!> module FILE I/O (netcdf)
!!
!! @par Description
!!          general file I/O module
!!          frontend interface of netCDF I/O routine
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-04-28 (R.Yoshida)   [new]
!!
!<
module mod_realinput
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  use scale_grid_real, only: &
     LON  => REAL_LON,  &
     LAT  => REAL_LAT,  &
     LONX => REAL_LONX, &
     LATY => REAL_LATY, &
     CZ   => REAL_CZ,   &
     FZ   => REAL_FZ
  use scale_index
  use scale_tracer
  use gtool_file_h
  use gtool_file, only: &
     FileGetShape, &
     FileRead
  use scale_process, only: &
     myrank => PRC_myrank,  &
     PRC_master,            &
     PRC_MPIstop
  use scale_const, only: &
     PI     => CONST_PI,     &
     Rdry   => CONST_Rdry,   &
     CPdry  => CONST_CPdry,  &
     CVdry  => CONST_CVdry,  &
     Rvap   => CONST_Rvap,   &
     CPvap  => CONST_CPvap,  &
     CVvap  => CONST_CVvap,  &
     CL     => CONST_CL,     &
     CI     => CONST_CI,     &
     PRE00  => CONST_PRE00,  &
     TEM00  => CONST_TEM00,  &
     grav   => CONST_GRAV,   &
     r_in_m => CONST_RADIUS, &
     dwatr  => CONST_DWATR,  &
     I_LW   => CONST_I_LW,   &
     I_SW   => CONST_I_SW
  use scale_external_io
  use scale_land_grid_index
  use scale_land_grid, only: &
     LCZ  => GRID_LCZ
  use scale_comm, only: &
     COMM_bcast
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtomSetup
  public :: ParentAtomInput
  public :: ParentAtomBoundary
  public :: ParentSurfaceInput

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: InputAtomSCALE
  private :: InputAtomWRF
  private :: InputAtomNICAM
  private :: InputSurfaceSCALE
  private :: InputSurfaceWRF
  private :: InputSurfaceNICAM
  private :: latlonz_interpolation_fact
  private :: interp_OceanLand_data
  private :: diagnose_number_concentration
  private :: haversine

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: iSCALE  = 1
  integer, private, parameter :: iWRFARW = 2
  integer, private, parameter :: iNICAM  = 3
  integer, private, parameter :: iJMAMSM = 4

  real(RP), private, parameter :: large_number_one   = 9.999E+15_RP
  real(RP), private, parameter :: large_number_two   = 8.888E+15_RP
  real(RP), private, parameter :: large_number_three = 7.777E+15_RP
  integer, parameter :: itp_nh = 3
  integer, parameter :: itp_nv = 2

  integer, private :: interp_search_divnum = 10
  logical, private :: wrfout = .false.  ! file type switch (wrfout or wrfrst)

  integer, parameter :: cosin = 1
  integer, parameter :: sine  = 2

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ParentAtomSetup( &
      dims,             &
      timelen,          &
      mdlid,            &
      basename_org,     &
      filetype,         &
      search_divnum_in, &
      wrf_file_type_in  )
    implicit none

    integer,          intent(out) :: dims(:)
    integer,          intent(out) :: timelen
    integer,          intent(out) :: mdlid
    character(len=*), intent(in)  :: basename_org
    character(len=*), intent(in)  :: filetype
    integer,          intent(in)  :: search_divnum_in
    logical,          intent(in)  :: wrf_file_type_in

    integer :: dims_ncm(4)
    character(len=H_LONG) :: basename
    intrinsic shape
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Setup]'

    select case(trim(filetype))
    case('SCALE-LES')
       if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: SCALE-LES'
       mdlid = iSCALE
       dims(:) = 0 ! unused
       timelen = 1 !temtative approach

    case('WRF-ARW')
       if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: WRF-ARW'
       mdlid = iWRFARW
       call ExternalFileGetShape( dims(:), timelen, mdlid, basename_org, myrank, single=.true. )
       if ( wrf_file_type_in ) then
          wrfout = .true.
          if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF History Output'
       else
          wrfout = .false.
          if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF Restart'
       endif

    case('NICAM-NETCDF')
       if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: NICAM-NETCDF'
       mdlid = iNICAM
       basename = "ms_pres"//trim(basename_org)
       call FileGetShape( dims_ncm(:), trim(basename), "ms_pres", 1, single=.true. )
       timelen = dims_ncm(4)
       dims(:) = 0
       dims(1) = dims_ncm(1)
       dims(2) = dims_ncm(2)
       dims(3) = dims_ncm(3)
       dims(4) = dims_ncm(1) ! nicam lat-lon data doesn't have staggered grid system
       dims(5) = dims_ncm(2) ! nicam lat-lon data doesn't have staggered grid system
       dims(6) = dims_ncm(3) ! nicam lat-lon data doesn't have staggered grid system
       basename = "la_tg"//trim(basename_org)
       call FileGetShape( dims_ncm(:), trim(basename), "la_tg", 1, single=.true. )
       dims(7) = dims_ncm(3) ! vertical grid of land model [tentative]

    !case('JMA-MSM')

    !case('NCEP-FNL')

    case default
       write(*,*) ' xxx Unsupported FILE TYPE:', trim(filetype)
       call PRC_MPIstop
    endselect

    interp_search_divnum = search_divnum_in
    if( IO_L ) write(IO_FID_LOG,*) '+++ Interpolation Search Block Dividing Num:', &
                                    interp_search_divnum

    return
  end subroutine ParentAtomSetup

  !-----------------------------------------------------------------------------
  !> Atmosphere Data Read
  subroutine ParentAtomInput( &
      dens,             &
      momz,             &
      momx,             &
      momy,             &
      rhot,             &
      qtrc,             &
      basename_org,     &
      use_file_density, &
      dims,             &
      mdlid,            &
      mptype_parent,    &
      timelen,          &
      serial            )
    implicit none

    real(RP),         intent(out) :: dens(:,:,:,:)
    real(RP),         intent(out) :: momz(:,:,:,:)
    real(RP),         intent(out) :: momx(:,:,:,:)
    real(RP),         intent(out) :: momy(:,:,:,:)
    real(RP),         intent(out) :: rhot(:,:,:,:)
    real(RP),         intent(out) :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in)  :: basename_org
    logical,          intent(in)  :: use_file_density ! use density data from files
    integer,          intent(in)  :: dims(:)
    integer,          intent(in)  :: mdlid            ! model type id
    integer,          intent(in)  :: mptype_parent    ! microphysics type of the parent model (number of classes)
    integer,          intent(in)  :: timelen          ! time steps in one file
    logical,          intent(in)  :: serial           ! read by a serial process

    integer :: k, i, j
    character(len=H_SHORT)  :: mptype_run
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Input]'

    select case (TRACER_TYPE)
    case ("DRY")
       mptype_run = 'dry'
    case ("KESSLER")
       mptype_run = 'single'
    case ("TOMITA08")
       mptype_run = 'single'
    case ("SN14")
       mptype_run = 'double'
    case ("SUZUKI10")
       mptype_run = 'double'
    case default
       write(*,*) 'xxx Unsupported TRACER_TYPE (', trim(TRACER_TYPE), '). Check!'
       call PRC_MPIstop
    end select

    if( mdlid == iSCALE ) then ! TYPE: SCALE-LES
       call InputAtomSCALE( dens    (:,:,:,:),   &
                            momz    (:,:,:,:),   &
                            momx    (:,:,:,:),   &
                            momy    (:,:,:,:),   &
                            rhot    (:,:,:,:),   &
                            qtrc    (:,:,:,:,:), &
                            basename_org,        &
                            use_file_density,    &
                            1,                   &  ! start time
                            timelen              )  ! end time

    elseif( mdlid == iWRFARW ) then ! TYPE: WRF-ARW
       call InputAtomWRF( dens    (:,:,:,:),   &
                          momz    (:,:,:,:),   &
                          momx    (:,:,:,:),   &
                          momy    (:,:,:,:),   &
                          rhot    (:,:,:,:),   &
                          qtrc    (:,:,:,:,:), &
                          basename_org,        &
                          dims(:),             &
                          mdlid,               &
                          mptype_run,          &
                          mptype_parent,       &
                          1,                   &  !suffix of start time: ts
                          timelen,             &  !suffix of end time: te
                          serial               )

    elseif( mdlid == iNICAM ) then ! TYPE: NICAM-NETCDF
       call InputAtomNICAM( dens    (:,:,:,:),   &
                            momz    (:,:,:,:),   &
                            momx    (:,:,:,:),   &
                            momy    (:,:,:,:),   &
                            rhot    (:,:,:,:),   &
                            qtrc    (:,:,:,:,:), &
                            basename_org,        &
                            dims(:),             &
                            1,                   &  ! start time
                            timelen              )  ! end time

    endif

    return
  end subroutine ParentAtomInput

  !-----------------------------------------------------------------------------
  !> Boundary Data Write
  subroutine ParentAtomBoundary( &
      dens,      &
      momz,      &
      momx,      &
      momy,      &
      rhot,      &
      qtrc,      &
      numsteps,  &
      update_dt, &
      basename,  &
      title      )
    use scale_comm, only: &
       COMM_vars, &
       COMM_wait
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    real(RP),         intent(in)   :: dens(:,:,:,:)
    real(RP),         intent(in)   :: momz(:,:,:,:)
    real(RP),         intent(in)   :: momx(:,:,:,:)
    real(RP),         intent(in)   :: momy(:,:,:,:)
    real(RP),         intent(in)   :: rhot(:,:,:,:)
    real(RP),         intent(in)   :: qtrc(:,:,:,:,:)
    real(RP),         intent(in)   :: update_dt
    character(LEN=*), intent(in)   :: basename
    character(LEN=*), intent(in)   :: title
    integer,          intent(in)   :: numsteps ! total time steps

    character(len=H_MID)  :: atmos_boundary_out_dtype = 'DEFAULT'  !< REAL4 or REAL8
    real(RP), allocatable :: buffer(:,:,:,:)

    integer :: k, i, j, n, iq
    integer :: ts, te
    integer :: timetarg
    !---------------------------------------------------------------------------

    ts = 1
    te = numsteps

    allocate( buffer(KA,IA,JA,te-ts+1) )

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Boundary]'

    call FILEIO_write( DENS(:,:,:,ts:te), basename, title,           &
                       'DENS', 'Reference Density', 'kg/m3', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt           )

    do n = ts, te
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMZ(k,i,j,n) / ( DENS(k+1,i,j,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer, basename, title,                 &
                       'VELZ', 'Reference VELZ', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt      )

    do n = ts, te
    do j = 1, JA
    do i = 1, IA-1
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMX(k,i,j,n) / ( DENS(k,i+1,j,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    do n = ts, te
       buffer(:,IA,:,n-ts+1) = buffer(:,IA-1,:,n-ts+1)
    end do
    call FILEIO_write( buffer, basename, title,                 &
                       'VELX', 'Reference VELX', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt      )

    do n = ts, te
    do j = 1, JA-1
    do i = 1, IA
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMY(k,i,j,n) / ( DENS(k,i,j+1,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    do n = ts, te
       buffer(:,:,JA,n-ts+1) = buffer(:,:,JA-1,n-ts+1)
    end do
    call FILEIO_write( buffer, basename, title,                 &
                       'VELY', 'Reference VELY', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt      )

    do n = ts, te
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = RHOT(k,i,j,n) / DENS(k,i,j,n)
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer, basename, title,                       &
                       'POTT', 'Reference PT', 'K', 'ZXYT',           &
                       atmos_boundary_out_dtype, update_dt            )

    do iq = 1, QA
       do n = ts, te
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          buffer(k,i,j,n-ts+1) = QTRC(k,i,j,n,iq)
       end do
       end do
       end do
       end do
       call FILEIO_write( buffer, basename, title,                                           &
                          AQ_NAME(iq), 'Reference '//trim(AQ_NAME(iq)), AQ_UNIT(iq), 'ZXYT', &
                          atmos_boundary_out_dtype, update_dt                                )
    end do

    deallocate( buffer )

    return
  end subroutine ParentAtomBoundary

  !-----------------------------------------------------------------------------
  !> Land&Ocean Data Read/Write
  !-----------------------------------------------------------------------------
  subroutine ParentSurfaceInput( &
      dens,                 &
      momz,                 &
      momx,                 &
      momy,                 &
      rhot,                 &
      qtrc,                 &
      basename_org,         &
      dims,                 &
      timelen,              &
      skiplen,              &
      use_file_landwater,   &
      init_landwater_ratio, &
      mdlid,                &
      serial                )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_ocean_vars, only: &
       OCEAN_vars_restart_write, &
       OCEAN_vars_external_in
    use mod_land_vars, only: &
       LAND_vars_restart_write, &
       LAND_vars_external_in
    use mod_urban_vars, only: &
       URBAN_vars_restart_write, &
       URBAN_vars_external_in
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_external_in
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    implicit none

    real(RP),         intent(in) :: dens(:,:,:,:)
    real(RP),         intent(in) :: momz(:,:,:,:)
    real(RP),         intent(in) :: momx(:,:,:,:)
    real(RP),         intent(in) :: momy(:,:,:,:)
    real(RP),         intent(in) :: rhot(:,:,:,:)
    real(RP),         intent(in) :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in) :: basename_org
    integer,          intent(in) :: dims(:)
    integer,          intent(in) :: timelen              ! time steps in one file
    integer,          intent(in) :: skiplen              ! skipped time steps
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    real(RP),         intent(in) :: init_landwater_ratio ! Ratio of land water to storage is constant,
                                                         !          if use_file_landwater is ".false."
    integer,          intent(in) :: mdlid                ! model type id
    logical,          intent(in) :: serial               ! read by a serial process

    real(RP) :: temp(KA,IA,JA)
    real(RP) :: pres(KA,IA,JA)

    real(RP) :: tg   (LKMAX,IA,JA)
    real(RP) :: strg (LKMAX,IA,JA)
    real(RP) :: roff (IA,JA)
    real(RP) :: qvef (IA,JA)
    real(RP) :: tw   (IA,JA)
    real(RP) :: lst  (IA,JA)
    real(RP) :: ust  (IA,JA)
    real(RP) :: sst  (IA,JA)
    real(RP) :: albw (IA,JA,2)  ! albedo for ocean
    real(RP) :: albg (IA,JA,2)  ! albedo for land
    real(RP) :: albu (IA,JA,2)  ! albedo for urban
    real(RP) :: z0w  (IA,JA)
    real(RP) :: skint(IA,JA)
    real(RP) :: skinw(IA,JA)
    real(RP) :: snowq(IA,JA)
    real(RP) :: snowt(IA,JA)

    real(RP) :: tc_urb(IA,JA)
    real(RP) :: qc_urb(IA,JA)
    real(RP) :: uc_urb(IA,JA)

    real(RP) :: init_qtrc(KA,IA,JA,QA)

    integer :: k, i, j, iq

    intrinsic shape
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Input-Surface]'

    if( LKMAX < 4 )then
       write(*,*) 'xxx LKMAX less than 4: ', LKMAX
       write(*,*) 'xxx in Real Case, LKMAX should be set more than 4'
       call PRC_MPIstop
    endif

    ! Read Data
    if( mdlid == iSCALE ) then ! TYPE: SCALE-LES
       call InputSurfaceSCALE( tg   (:,:,:),         &
                               strg (:,:,:),         &
                               roff (:,:),           &
                               qvef (:,:),           &
                               tw   (:,:),           &
                               lst  (:,:),           &
                               ust  (:,:),           &
                               sst  (:,:),           &
                               albw (:,:,:),         &
                               albg (:,:,:),         &
                               z0w  (:,:),           &
                               skint(:,:),           &
                               skinw(:,:),           &
                               snowq(:,:),           &
                               snowt(:,:),           &
                               basename_org,         &
                               use_file_landwater,   &  ! use land water data from files
                               init_landwater_ratio  )  ! Ratio of land water to storage is constant
    elseif( mdlid == iWRFARW ) then ! TYPE: WRF-ARW
       call InputSurfaceWRF( tg   (:,:,:),         &
                             strg (:,:,:),         &
                             roff (:,:),           &
                             qvef (:,:),           &
                             tw   (:,:),           &
                             lst  (:,:),           &
                             ust  (:,:),           &
                             sst  (:,:),           &
                             albw (:,:,:),         &
                             albg (:,:,:),         &
                             z0w  (:,:),           &
                             skint(:,:),           &
                             skinw(:,:),           &
                             snowq(:,:),           &
                             snowt(:,:),           &
                             basename_org,         &
                             dims(:),              &
                             use_file_landwater,   &  ! use land water data from files
                             init_landwater_ratio, &  ! Ratio of land water to storage is constant
                             mdlid,                &
                             serial                )
    elseif( mdlid == iNICAM ) then ! TYPE: NICAM-NETCDF
       call InputSurfaceNICAM( tg   (:,:,:),        &
                               strg (:,:,:),        &
                               roff (:,:),          &
                               qvef (:,:),          &
                               tw   (:,:),          &
                               lst  (:,:),          &
                               ust  (:,:),          &
                               sst  (:,:),          &
                               albw (:,:,:),        &
                               albg (:,:,:),        &
                               z0w  (:,:),          &
                               skint(:,:),          &
                               skinw(:,:),          &
                               snowq(:,:),          &
                               snowt(:,:),          &
                               basename_org,        &
                               dims(:),             &
                               skiplen,             &   ! the number of skipped data
                               use_file_landwater,  &   ! use land water data from files
                               init_landwater_ratio )   ! Ratio of land water to storage is constant
    endif

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
    do iq = 1, QA
       init_qtrc(k,i,j,iq) = qtrc(k,i,j,1,iq)
    end do
    end do
    end do
    end do

    call THERMODYN_temp_pres( temp     (:,:,:),   & ! [OUT]
                              pres     (:,:,:),   & ! [OUT] not used
                              dens     (:,:,:,1), & ! [IN]
                              rhot     (:,:,:,1), & ! [IN]
                              init_qtrc(:,:,:,:)  ) ! [IN]

    do j = 1, JA
    do i = 1, IA
       tc_urb(i,j) = temp(KS,i,j)
#ifdef DRY
       qc_urb(i,j) = 0.0_RP
#else
       qc_urb(i,j) = qtrc(KS,i,j,1,I_QV)
#endif
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA-1
       uc_urb(i,j) = max(sqrt( ( momx(KS,i,j,1) / ( dens(KS,i+1,  j,1) + dens(KS,i,j,1) ) * 2.0_RP )**2.0_RP &
                             + ( momy(KS,i,j,1) / ( dens(KS,  i,j+1,1) + dens(KS,i,j,1) ) * 2.0_RP )**2.0_RP ), 0.01_RP)
    enddo
    enddo
    do j = 1, JA-1
       uc_urb(IA,j) = max(sqrt( ( momx(KS,IA,j,1) / dens(KS,IA,j,1) )**2.0_RP &
                              + ( momy(KS,IA,j,1) / ( dens(KS,IA,j+1,1) + dens(KS,IA,j,1) ) * 2.0_RP )**2.0_RP ), 0.01_RP)
    enddo
    do i = 1, IA-1
       uc_urb(i,JA) = max(sqrt( ( momx(KS,i,JA,1) / ( dens(KS,i+1,JA,1) + dens(KS,i,JA,1) ) * 2.0_RP )**2.0_RP &
                              + ( momy(KS,i,JA,1) / dens(KS,i,JA,1) )**2.0_RP ), 0.01_RP)
    enddo
    uc_urb(IA,JA) = max(sqrt( ( momx(KS,IA,JA,1) / dens(KS,IA,JA,1) )**2.0_RP &
                            + ( momy(KS,IA,JA,1) / dens(KS,IA,JA,1) )**2.0_RP ), 0.01_RP)

    ! copy albedo of land to urban
    do j = 1, JA
    do i = 1, IA
       albu(i,j,:) = albg(i,j,:)
    enddo
    enddo

    call COMM_vars8( uc_urb, 1 )
    call COMM_wait ( uc_urb, 1, .false. )

    ! Write out: Ocean
    call OCEAN_vars_external_in( tw, sst, albw, z0w, z0w, z0w )

    ! Write out: Land
    call LAND_vars_external_in( tg, strg, lst, albg )

    ! Write out: Urban
    call URBAN_vars_external_in( tc_urb, qc_urb, uc_urb, ust, albu )

    ! Input to PHY_SF Container
    call ATMOS_PHY_SF_vars_external_in( skint, albw(:,:,I_LW), albw(:,:,I_SW), z0w, z0w, z0w )

    return
  end subroutine ParentSurfaceInput

  !-----------------------------------------------------------------------------
  !> Individual Procedures
  !-----------------------------------------------------------------------------
  subroutine InputAtomSCALE( &
      dens,             & ! (out)
      momz,             & ! (out)
      momx,             & ! (out)
      momy,             & ! (out)
      rhot,             & ! (out)
      qtrc,             & ! (out)
      basename_org,     & ! (in)
      use_file_density, & ! (in)
      start_step,       & ! (in)
      end_step          ) ! (in)
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_grid_nest, only: &
       PARENT_KMAX,     &
       PARENT_IMAX,     &
       PARENT_JMAX,     &
       NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y, &
       NEST_TILE_ID
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
       THERMODYN_pott      => ATMOS_THERMODYN_pott
    use scale_gridtrans, only: &
       rotc => GTRANS_ROTC
    use scale_topography, only: &
       topo => TOPO_Zsfc
    implicit none

    real(RP),         intent(out) :: dens(:,:,:,:)
    real(RP),         intent(out) :: momz(:,:,:,:)
    real(RP),         intent(out) :: momx(:,:,:,:)
    real(RP),         intent(out) :: momy(:,:,:,:)
    real(RP),         intent(out) :: rhot(:,:,:,:)
    real(RP),         intent(out) :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in)  :: basename_org
    logical,          intent(in)  :: use_file_density
    integer,          intent(in)  :: start_step
    integer,          intent(in)  :: end_step

    integer, parameter    :: handle = 1

    ! work
    real(RP), allocatable :: read2D(:,:)
    real(RP), allocatable :: read3D(:,:,:)

    real(RP), allocatable :: lon_org (:,:,:)
    real(RP), allocatable :: lat_org (:,:,:)
    real(RP), allocatable :: lonu_org(:,:,:)
    real(RP), allocatable :: latu_org(:,:,:)
    real(RP), allocatable :: lonv_org(:,:,:)
    real(RP), allocatable :: latv_org(:,:,:)
    real(RP), allocatable :: geoh_org(:,:,:,:)
    real(RP), allocatable :: geof_org(:,:,:,:)

    real(RP), allocatable :: dens_org(:,:,:,:)
    real(RP), allocatable :: momz_org(:,:,:,:)
    real(RP), allocatable :: momx_org(:,:,:,:)
    real(RP), allocatable :: momy_org(:,:,:,:)
    real(RP), allocatable :: rhot_org(:,:,:,:)
    real(RP), allocatable :: qtrc_org(:,:,:,:,:)

    real(RP), allocatable :: tsfc_org(:,:,:)
    real(RP), allocatable :: qsfc_org(:,:,:,:)
    real(RP), allocatable :: mslp_org(:,:,:)

    real(RP), allocatable :: velz_org(:,:,:,:)
    real(RP), allocatable :: velx_org(:,:,:,:)
    real(RP), allocatable :: vely_org(:,:,:,:)
    real(RP), allocatable :: pott_org(:,:,:,:)
    real(RP), allocatable :: temp_org(:,:,:,:)
    real(RP), allocatable :: pres_org(:,:,:,:)

    real(RP) :: velz  (KA,IA,JA,start_step:end_step)
    real(RP) :: velx  (KA,IA,JA,start_step:end_step)
    real(RP) :: vely  (KA,IA,JA,start_step:end_step)
    real(RP) :: llvelx(KA,IA,JA,start_step:end_step)
    real(RP) :: llvely(KA,IA,JA,start_step:end_step)
    real(RP) :: work  (KA,IA,JA,start_step:end_step)
    real(RP) :: pott  (KA,IA,JA,start_step:end_step)
    real(RP) :: temp  (KA,IA,JA,start_step:end_step)
    real(RP) :: pres  (KA,IA,JA,start_step:end_step)

    real(RP) :: pott_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: pres_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: temp_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: qtrc_sfc(1,IA,JA,start_step:end_step,QA)
    real(RP) :: mslp_sfc(1,IA,JA,start_step:end_step)

    real(RP) :: hfact(   IA,JA,itp_nh       )
    real(RP) :: vfact(KA,IA,JA,itp_nh,itp_nv)

    integer :: igrd(   IA,JA,itp_nh       )
    integer :: jgrd(   IA,JA,itp_nh       )
    integer :: kgrd(KA,IA,JA,itp_nh,itp_nv)

    real(RP) :: qc(KA,IA,JA)
    real(RP) :: qc_sfc(1,IA,JA)

    real(RP) :: wgt_up, wgt_bt
    real(RP) :: z1, z2
    real(RP) :: pres1, pres2

    integer :: KALL, IALL, JALL
    integer :: rank

    integer :: k, i, j, n, iq
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye

    logical :: lack_of_val
    !---------------------------------------------------------------------------

    KALL = PARENT_KMAX(handle)
    IALL = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    JALL = PARENT_JMAX(handle) * NEST_TILE_NUM_Y

    allocate( read2D( PARENT_IMAX(handle), PARENT_JMAX(handle)                      ) )
    allocate( read3D( PARENT_IMAX(handle), PARENT_JMAX(handle), PARENT_KMAX(handle) ) )

    allocate( lon_org (       IALL, JALL, start_step:end_step )    )
    allocate( lat_org (       IALL, JALL, start_step:end_step )    )
    allocate( lonu_org(       IALL, JALL, start_step:end_step )    )
    allocate( latu_org(       IALL, JALL, start_step:end_step )    )
    allocate( lonv_org(       IALL, JALL, start_step:end_step )    )
    allocate( latv_org(       IALL, JALL, start_step:end_step )    )
    allocate( geoh_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( geof_org( KALL, IALL, JALL, start_step:end_step )    )

    allocate( dens_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( momz_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( momx_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( momy_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( rhot_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( qtrc_org( KALL, IALL, JALL, start_step:end_step, QA) )

    allocate( tsfc_org(       IALL, JALL, start_step:end_step )    )
    allocate( qsfc_org(       IALL, JALL, start_step:end_step, QA) )
    allocate( mslp_org(       IALL, JALL, start_step:end_step )    )

    allocate( velz_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( velx_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( vely_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( pott_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( temp_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( pres_org( KALL, IALL, JALL, start_step:end_step )    )

    if( IO_L ) write(IO_FID_LOG,*) ''
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

       call FileRead( read2D(:,:),   BASENAME_ORG, "lon",        1, rank )
       lon_org (xs:xe,ys:ye,1)  = read2D(:,:) * D2R

       call FileRead( read2D(:,:),   BASENAME_ORG, "lat",        1, rank )
       lat_org (xs:xe,ys:ye,1)  = read2D(:,:) * D2R

       call FileRead( read2D(:,:),   BASENAME_ORG, "lon_uy",     1, rank )
       lonu_org (xs:xe,ys:ye,1) = read2D(:,:) * D2R

       call FileRead( read2D(:,:),   BASENAME_ORG, "lat_uy",     1, rank )
       latu_org (xs:xe,ys:ye,1) = read2D(:,:) * D2R

       call FileRead( read2D(:,:),   BASENAME_ORG, "lon_xv",     1, rank )
       lonv_org (xs:xe,ys:ye,1) = read2D(:,:) * D2R

       call FileRead( read2D(:,:),   BASENAME_ORG, "lat_xv",     1, rank )
       latv_org (xs:xe,ys:ye,1) = read2D(:,:) * D2R

       call FileRead( read3D(:,:,:), BASENAME_ORG, "height",     1, rank )
       do k = 1, KALL
         geoh_org(k,xs:xe,ys:ye,1) = read3D(:,:,k)
       end do

       call FileRead( read3D(:,:,:), BASENAME_ORG, "height_xyw", 1, rank )
       do k = 1, KALL
         geof_org(k,xs:xe,ys:ye,1) = read3D(:,:,k)
       end do

       do n = start_step, end_step
         call FileRead( read2D(:,:), BASENAME_ORG, "T2", n, rank )
         tsfc_org(xs:xe,ys:ye,n) = read2D(:,:)
       end do

       do n = start_step, end_step
         call FileRead( read2D(:,:), BASENAME_ORG, "Q2", n, rank )
         qsfc_org(xs:xe,ys:ye,n,I_QV) = read2D(:,:)
       end do

       do n = start_step, end_step
         call FileRead( read2D(:,:), BASENAME_ORG, "MSLP", n, rank )
         mslp_org(xs:xe,ys:ye,n) = read2D(:,:)
       end do

       do n = start_step, end_step
         call FileRead( read3D(:,:,:), BASENAME_ORG, "DENS", n, rank )
         do k = 1, KALL
           dens_org(k,xs:xe,ys:ye,n) = read3D(:,:,k)
         end do
       end do

       do n = start_step, end_step
         call FileRead( read3D(:,:,:), BASENAME_ORG, "MOMZ", n, rank )
         do k = 1, KALL
           momz_org(k,xs:xe,ys:ye,n) = read3D(:,:,k)
         end do
       end do

       do n = start_step, end_step
         call FileRead( read3D(:,:,:), BASENAME_ORG, "MOMX", n, rank )
         do k = 1, KALL
           momx_org(k,xs:xe,ys:ye,n) = read3D(:,:,k)
         end do
       end do

       do n = start_step, end_step
         call FileRead( read3D(:,:,:), BASENAME_ORG, "MOMY", n, rank )
         do k = 1, KALL
           momy_org(k,xs:xe,ys:ye,n) = read3D(:,:,k)
         end do
       end do

       do n = start_step, end_step
         call FileRead( read3D(:,:,:), BASENAME_ORG, "RHOT", n, rank )
         do k = 1, KALL
           rhot_org(k,xs:xe,ys:ye,n) = read3D(:,:,k)
         end do
       end do

       do iq = 1, QA
         do n = start_step, end_step
           call FileRead( read3D(:,:,:), BASENAME_ORG, AQ_NAME(iq), n, rank )
           do k = 1, KALL
             qtrc_org(k,xs:xe,ys:ye,n,iq) = read3D(:,:,k)
           end do
         end do
       end do

    end do

    ! fill all step
    do n = start_step, end_step
      lon_org (:,:,n)   = lon_org (:,:,1)
      lat_org (:,:,n)   = lat_org (:,:,1)
      lonu_org(:,:,n)   = lonu_org(:,:,1)
      latu_org(:,:,n)   = latu_org(:,:,1)
      lonv_org(:,:,n)   = lonv_org(:,:,1)
      latv_org(:,:,n)   = latv_org(:,:,1)
      geoh_org(:,:,:,n) = geoh_org(:,:,:,1)
      geof_org(:,:,:,n) = geof_org(:,:,:,1)
    end do

    call check_domain_compatibility( lon_org(:,:,1),  lat_org(:,:,1),  geoh_org(:,:,:,1), &
                                     LON(:,:),        LAT(:,:),        CZ(:,:,:) )
    call check_domain_compatibility( lonu_org(:,:,1), latu_org(:,:,1), geoh_org(:,:,:,1), &
                                     LONX(:,:),       LAT(:,:),        CZ(:,:,:), skip_z=.true. )
    call check_domain_compatibility( lonv_org(:,:,1), latv_org(:,:,1), geoh_org(:,:,:,1), &
                                     LON(:,:),        LATY(:,:),       CZ(:,:,:), skip_z=.true. )
    call check_domain_compatibility( lon_org(:,:,1),  lat_org(:,:,1),  geof_org(:,:,:,1), &
                                     LON(:,:),        LAT(:,:),        FZ(:,:,:), skip_x=.true., skip_y=.true. )

    ! for vector (w) points
    call latlonz_interpolation_fact( hfact   (:,:,:),     & ! [OUT]
                                     vfact   (:,:,:,:,:), & ! [OUT]
                                     kgrd    (:,:,:,:,:), & ! [OUT]
                                     igrd    (:,:,:),     & ! [OUT]
                                     jgrd    (:,:,:),     & ! [OUT]
                                     FZ      (:,:,:),     & ! [IN]
                                     LAT     (:,:),       & ! [IN]
                                     LON     (:,:),       & ! [IN]
                                     geof_org(:,:,:,:),   & ! [IN]
                                     lat_org (:,:,:),     & ! [IN]
                                     lon_org (:,:,:),     & ! [IN]
                                     KALL, IALL, JALL,    & ! [IN]
                                     start_step           ) ! [IN]

    ! convert from momentum to velocity
    do n = start_step, end_step
    do j = 1, JALL
    do i = 1, IALL
    do k = 1, KALL-1
       velz_org(k,i,j,n) = momz_org(k,i,j,n) / ( dens_org(k+1,i,j,n) + dens_org(k,i,j,n) ) * 2.0_RP
    end do
    end do
    end do
    end do
    do n = start_step, end_step
    do j = 1, JALL
    do i = 1, IALL
       velz_org(KALL,i,j,n) = 0.0_RP
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE-1
       velz(k,i,j,n) = velz_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + velz_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + velz_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + velz_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + velz_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + velz_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    end do
    end do
    end do
    end do

    ! for vector (u) points
    call latlonz_interpolation_fact( hfact   (:,:,:),     & ! [OUT]
                                     vfact   (:,:,:,:,:), & ! [OUT]
                                     kgrd    (:,:,:,:,:), & ! [OUT]
                                     igrd    (:,:,:),     & ! [OUT]
                                     jgrd    (:,:,:),     & ! [OUT]
                                     CZ      (:,:,:),     & ! [IN]
                                     LAT     (:,:),       & ! [IN]
                                     LONX    (:,:),       & ! [IN]
                                     geoh_org(:,:,:,:),   & ! [IN]
                                     latu_org(:,:,:),     & ! [IN]
                                     lonu_org(:,:,:),     & ! [IN]
                                     KALL, IALL, JALL,    & ! [IN]
                                     start_step           ) ! [IN]

    ! convert from momentum to velocity
    do n = start_step, end_step
    do j = 1, JALL
    do i = 1, IALL-1
    do k = 1, KALL
       velx_org(k,i,j,n) = momx_org(k,i,j,n) / ( dens_org(k,i+1,j,n) + dens_org(k,i,j,n) ) * 2.0_RP
    end do
    end do
    end do
    end do
    do n = start_step, end_step
    do j = 1, JALL
    do k = 1, KALL
       velx_org(k,IALL,j,n) = momx_org(k,IALL,j,n) / dens_org(k,IALL,j,n)
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
       work(k,i,j,n) = velx_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + velx_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + velx_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + velx_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + velx_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + velx_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    end do
    end do
    end do
    end do
    ! from staggered point to scalar point
    do n = start_step, end_step
       do j = 1, JA
       do i = 2, IA
       do k = KS-1, KE+1
          llvelx(k,i,j,n) = ( work(k,i-1,j,n) + work(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do k = KS-1, KE+1
          llvelx(k,1,j,n) = work(k,1,j,n)
       end do
       end do
       call COMM_vars8( llvelx(:,:,:,n), 1 )
       call COMM_wait ( llvelx(:,:,:,n), 1, .false. )
    end do

    ! for vector (v) points
    call latlonz_interpolation_fact( hfact   (:,:,:),     & ! [OUT]
                                     vfact   (:,:,:,:,:), & ! [OUT]
                                     kgrd    (:,:,:,:,:), & ! [OUT]
                                     igrd    (:,:,:),     & ! [OUT]
                                     jgrd    (:,:,:),     & ! [OUT]
                                     CZ      (:,:,:),     & ! [IN]
                                     LATY    (:,:),       & ! [IN]
                                     LON     (:,:),       & ! [IN]
                                     geoh_org(:,:,:,:),   & ! [IN]
                                     latv_org(:,:,:),     & ! [IN]
                                     lonv_org(:,:,:),     & ! [IN]
                                     KALL, IALL, JALL,    & ! [IN]
                                     start_step           ) ! [IN]

    ! convert from momentum to velocity
    do n = start_step, end_step
    do j = 1, JALL-1
    do i = 1, IALL
    do k = 1, KALL
       vely_org(k,i,j,n) = momy_org(k,i,j,n) / ( dens_org(k,i,j+1,n) + dens_org(k,i,j,n) ) * 2.0_RP
    end do
    end do
    end do
    end do
    do n = start_step, end_step
    do i = 1, IALL
    do k = 1, KALL
       vely_org(k,i,JALL,n) = momy_org(k,i,JALL,n) / dens_org(k,i,JALL,n)
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
       work(k,i,j,n) = vely_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + vely_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + vely_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + vely_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + vely_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + vely_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    end do
    end do
    end do
    end do
    ! from staggered point to scalar point
    do n = start_step, end_step
       do j = 2, JA
       do i = 1, IA
       do k = KS-1, KE+1
          llvely(k,i,j,n) = ( work(k,i,j-1,n) + work(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS-1, KE+1
          llvely(k,i,1,n) = work(k,i,1,n)
       end do
       end do
       call COMM_vars8( llvely(:,:,:,n), 1 )
       call COMM_wait ( llvely(:,:,:,n), 1, .false. )
    end do

    do n = start_step, end_step
       ! convert from latlon coordinate to local mapping (x)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          work(k,i,j,n) = llvelx(k,i,j,n) * rotc(i,j,cosin) + llvely(k,i,j,n) * rotc(i,j,sine )
       end do
       end do
       end do

       ! from scalar point to staggered point
       do j = 1, JA
       do i = 1, IA-1
       do k = KS, KE
          velx(k,i,j,n) = ( work(k,i+1,j,n) + work(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do k = KS, KE
          velx(k,IA,j,n) = work(k,IA,j,n)
       end do
       end do
       velx(KS-1,:,:,n) = 0.0_RP
       velx(KS-2,:,:,n) = 0.0_RP
       call COMM_vars8( velx(:,:,:,n), n )
       call COMM_wait ( velx(:,:,:,n), n, .false. )

    end do

    do n = start_step, end_step
       ! convert from latlon coordinate to local mapping (y)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          work(k,i,j,n) = - llvelx(k,i,j,n) * rotc(i,j,sine ) + llvely(k,i,j,n) * rotc(i,j,cosin)
       end do
       end do
       end do

       do j = 1, JA-1
       do i = 1, IA
       do k = KS, KE
          vely(k,i,j,n) = ( work(k,i,j+1,n) + work(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          vely(k,i,JA,n) = work(k,i,JA,n)
       end do
       end do
       vely(KS-1,:,:,n) = 0.0_RP
       vely(KS-2,:,:,n) = 0.0_RP
       call COMM_vars8( vely(:,:,:,n), 1 )
       call COMM_wait ( vely(:,:,:,n), 1, .false. )

    end do

    ! for scalar points
    call latlonz_interpolation_fact( hfact   (:,:,:),     & ! [OUT]
                                     vfact   (:,:,:,:,:), & ! [OUT]
                                     kgrd    (:,:,:,:,:), & ! [OUT]
                                     igrd    (:,:,:),     & ! [OUT]
                                     jgrd    (:,:,:),     & ! [OUT]
                                     CZ      (:,:,:),     & ! [IN]
                                     LAT     (:,:),       & ! [IN]
                                     LON     (:,:),       & ! [IN]
                                     geoh_org(:,:,:,:),   & ! [IN]
                                     lat_org (:,:,:),     & ! [IN]
                                     lon_org (:,:,:),     & ! [IN]
                                     KALL, IALL, JALL,    & ! [IN]
                                     start_step           ) ! [IN]

    do n = start_step, end_step
    do j = 1, JALL
    do i = 1, IALL
    do k = 1, KALL
       do iq = 1, QA
          qtrc_org(k,i,j,n,iq) = max( qtrc_org(k,i,j,n,iq), 0.0_RP )
       end do

       ! diagnose temp and pres
       call THERMODYN_temp_pres( temp_org(k,i,j,n),   & ! [OUT]
                                 pres_org(k,i,j,n),   & ! [OUT]
                                 dens_org(k,i,j,n),   & ! [IN]
                                 rhot_org(k,i,j,n),   & ! [IN]
                                 qtrc_org(k,i,j,n,:)  ) ! [IN]

       pott_org(k,i,j,n) = rhot_org(k,i,j,n) / dens_org(k,i,j,n)
       dens_org(k,i,j,n) = log( dens_org(k,i,j,n) )
    end do
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
       pott(k,i,j,n) = pott_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + pott_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + pott_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + pott_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + pott_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + pott_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)

       pres(k,i,j,n) = pres_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + pres_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + pres_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + pres_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + pres_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + pres_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)

       do iq = 1, QA
          qtrc(k,i,j,n,iq) = qtrc_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n,iq) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                           + qtrc_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n,iq) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                           + qtrc_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n,iq) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                           + qtrc_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n,iq) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                           + qtrc_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n,iq) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                           + qtrc_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n,iq) * hfact(i,j,3) * vfact(k,i,j,3,2)
       end do
    end do
    end do
    end do
    end do

    if( use_file_density ) then
       ! use logarithmic density to interpolate more accurately
       do n = start_step, end_step
       do j = 1, JA
       do i = 1, IA
       do k = KS-1, KE+1
          dens(k,i,j,n) = exp( dens_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                             + dens_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                             + dens_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                             + dens_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                             + dens_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                             + dens_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2) )
       end do
       end do
       end do
       end do
    else
       do n = start_step, end_step
       do j = 1, JA
       do i = 1, IA
          temp_sfc(1,i,j,n) = tsfc_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                            + tsfc_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                            + tsfc_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)

          mslp_sfc(1,i,j,n) = mslp_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                            + mslp_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                            + mslp_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)

          ! Interpolate Surface pressure from SLP and PRES
          lack_of_val = .true.

          do k = KS-1, KE
             if( k == KS-1 ) then
               z1    = 0.0_RP
               z2    = CZ      (k+1,i,j  )
               pres1 = mslp_sfc(  1,i,j,n)
               pres2 = pres    (k+1,i,j,n)
             else
               z1    = CZ  (k  ,i,j  )
               z2    = CZ  (k+1,i,j  )
               pres1 = pres(k  ,i,j,n)
               pres2 = pres(k+1,i,j,n)
             endif
             if( topo(i,j) >= z1 .and. &
                 topo(i,j) <  z2       ) then
                lack_of_val = .false. ! found

                wgt_bt = ( z2        - topo(i,j) ) / (z2 - z1)
                wgt_up = ( topo(i,j) - z1        ) / (z2 - z1)

                pres_sfc(1,i,j,n) = exp( log(pres1)*wgt_bt + log(pres2)*wgt_up )

                exit ! exit loop
             endif
          enddo

          if( lack_of_val ) then
             write(IO_FID_LOG,*) 'realinput ATM SCALE: cannot estimate pres_sfc',i,j,n
             call PRC_MPIstop
          endif

          ! interpolate QV (=Q2) only: other QTRC are set zero
          qtrc_sfc(1,i,j,n,:)    = 0.0_RP
          qtrc_sfc(1,i,j,n,I_QV) = qsfc_org(igrd(i,j,1),jgrd(i,j,1),n,I_QV) * hfact(i,j,1) &
                                 + qsfc_org(igrd(i,j,2),jgrd(i,j,2),n,I_QV) * hfact(i,j,2) &
                                 + qsfc_org(igrd(i,j,3),jgrd(i,j,3),n,I_QV) * hfact(i,j,3)

          call THERMODYN_pott( pott_sfc(1,i,j,n),  & ! [OUT]
                               temp_sfc(1,i,j,n),  & ! [IN]
                               pres_sfc(1,i,j,n),  & ! [IN]
                               qtrc_sfc(1,i,j,n,:) ) ! [IN]
       end do
       end do
       end do

       do n = start_step, end_step
#ifdef DRY
          qc = 0.0_RP
          qc_sfc    = 0.0_RP
#else
          qc = qtrc(:,:,:,n,I_QC)
          qc_sfc = qtrc_sfc(:,:,:,n,I_QC)
#endif
          ! make density in moist condition
          call HYDROSTATIC_buildrho_real( dens    (:,:,:,n),      & ! [OUT]
                                          temp    (:,:,:,n),      & ! [OUT]
                                          pres    (:,:,:,n),      & ! [OUT]
                                          pott    (:,:,:,n),      & ! [IN]
                                          qtrc    (:,:,:,n,I_QV), & ! [IN]
                                          qc      (:,:,:),        & ! [OUT]
                                          temp_sfc(:,:,:,n),      & ! [OUT]
                                          pres_sfc(:,:,:,n),      & ! [IN]
                                          pott_sfc(:,:,:,n),      & ! [IN]
                                          qtrc_sfc(:,:,:,n,I_QV), & ! [IN]
                                          qc_sfc  (:,:,:)         ) ! [IN]

          call COMM_vars8( dens(:,:,:,n), 1 )
          call COMM_wait ( dens(:,:,:,n), 1 )
       end do
    end if

    do n = start_step, end_step
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE-1
          momz(k,i,j,n) = velz(k,i,j,n) * ( dens(k+1,i  ,j  ,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do i = 1, IA
          momz(KE,i,j,n) = 0.0_RP
       end do
       end do
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          rhot(k,i,j,n) = pott(k,i,j,n) * dens(k,i,j,n)
       end do
       end do
       end do
       do j = 1, JA
       do i = 1, IA-1
       do k = KS, KE
          momx(k,i,j,n) = velx(k,i,j,n) * ( dens(k  ,i+1,j  ,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do k = KS, KE
          momx(k,IA,j,n) = velx(k,IA,j,n) * dens(k,IA,j,n)
       end do
       end do
       call COMM_vars8( momx(:,:,:,n), 1 )

       do j = 1, JA-1
       do i = 1, IA
       do k = KS, KE
          momy(k,i,j,n) = vely(k,i,j,n) * ( dens(k  ,i  ,j+1,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          momy(k,i,JA,n) = vely(k,i,JA,n) * dens(k,i,JA,n)
       end do
       end do
       call COMM_vars8( momy(:,:,:,n), 2 )

       call COMM_wait ( momx(:,:,:,n), 1, .false. )
       call COMM_wait ( momy(:,:,:,n), 2, .false. )

    end do

    deallocate( read2D )
    deallocate( read3D )

    deallocate( lon_org  )
    deallocate( lat_org  )
    deallocate( lonu_org )
    deallocate( latu_org )
    deallocate( lonv_org )
    deallocate( latv_org )
    deallocate( geoh_org )
    deallocate( geof_org )

    deallocate( dens_org )
    deallocate( momz_org )
    deallocate( momx_org )
    deallocate( momy_org )
    deallocate( rhot_org )
    deallocate( qtrc_org )

    return
  end subroutine InputAtomSCALE

  !-----------------------------------------------------------------------------
  subroutine InputAtomWRF( &
      dens,          & ! (out)
      momz,          & ! (out)
      momx,          & ! (out)
      momy,          & ! (out)
      rhot,          & ! (out)
      qtrc,          & ! (out)
      basename,      & ! (in)
      dims,          & ! (in)
      mdlid,         & ! (in)
      mptype_run,    & ! (in)
      mptype_parent, & ! (in)
      ts,            & ! (in)
      te,            & ! (in)
      serial         ) ! (in)
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_atmos_thermodyn, only: &
       THERMODYN_pott => ATMOS_THERMODYN_pott
    use scale_gridtrans, only: &
       rotc => GTRANS_ROTC
    implicit none

    real(RP),         intent(out)  :: dens(:,:,:,:)
    real(RP),         intent(out)  :: momz(:,:,:,:)
    real(RP),         intent(out)  :: momx(:,:,:,:)
    real(RP),         intent(out)  :: momy(:,:,:,:)
    real(RP),         intent(out)  :: rhot(:,:,:,:)
    real(RP),         intent(out)  :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in)   :: basename
    character(LEN=*), intent(in)   :: mptype_run
    integer,          intent(in)   :: mptype_parent
    integer,          intent(in)   :: mdlid
    integer,          intent(in)   :: dims(:)
    integer,          intent(in)   :: ts      ! suffix of start time
    integer,          intent(in)   :: te      ! suffix of end time
    logical,          intent(in)   :: serial  ! read by a serial process

    ! Defined parameters in WRF
    real(RP), parameter :: t0  = 300.0_RP
    real(RP), parameter :: p0  = 1000.0E+2_RP
    real(RP), parameter :: Rd  = 287.04_RP
    real(RP), parameter :: Cp  = 7.0_RP * Rd / 2.0_RP
    real(RP), parameter :: RCP = Rd / Cp

    real(SP), allocatable :: read_xy (:,:,:)
    real(SP), allocatable :: read_uy (:,:,:)
    real(SP), allocatable :: read_xv (:,:,:)
    real(SP), allocatable :: read_zxy(:,:,:,:)
    real(SP), allocatable :: read_wxy(:,:,:,:)
    real(SP), allocatable :: read_zuy(:,:,:,:)
    real(SP), allocatable :: read_zxv(:,:,:,:)

    real(RP), allocatable :: velz_org(:,:,:,:)
    real(RP), allocatable :: velx_org(:,:,:,:)
    real(RP), allocatable :: vely_org(:,:,:,:)
    real(RP), allocatable :: pott_org(:,:,:,:)
    real(RP), allocatable :: temp_org(:,:,:,:)
    real(RP), allocatable :: pres_org(:,:,:,:)
    real(RP), allocatable :: p_org   (:,:,:,:)
    real(RP), allocatable :: pb_org  (:,:,:,:)
    real(RP), allocatable :: qtrc_org(:,:,:,:,:)

    real(RP), allocatable :: tsfc_org(:,:,:)
    real(RP), allocatable :: psfc_org(:,:,:)
    real(RP), allocatable :: qsfc_org(:,:,:,:)

    real(RP), allocatable :: ph_org  (:,:,:,:)
    real(RP), allocatable :: phb_org (:,:,:,:)
    real(RP), allocatable :: geof_org(:,:,:,:)
    real(RP), allocatable :: geoh_org(:,:,:,:)
    real(RP), allocatable :: lat_org (:,:,:)
    real(RP), allocatable :: lon_org (:,:,:)
    real(RP), allocatable :: latu_org(:,:,:)
    real(RP), allocatable :: lonu_org(:,:,:)
    real(RP), allocatable :: latv_org(:,:,:)
    real(RP), allocatable :: lonv_org(:,:,:)

    real(RP) :: velz  (KA,IA,JA)
    real(RP) :: velx  (KA,IA,JA)
    real(RP) :: vely  (KA,IA,JA)
    real(RP) :: llvelx(KA,IA,JA)
    real(RP) :: llvely(KA,IA,JA)
    real(RP) :: pott  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)
    real(RP) :: work  (KA,IA,JA)

    real(RP) :: pott_sfc(1,IA,JA)
    real(RP) :: temp_sfc(1,IA,JA)
    real(RP) :: pres_sfc(1,IA,JA)
    real(RP) :: qtrc_sfc(1,IA,JA,QA)

    real(RP) :: hfact(   IA,JA,itp_nh       )
    real(RP) :: vfact(KA,IA,JA,itp_nh,itp_nv)

    integer :: igrd(   IA,JA,itp_nh       )
    integer :: jgrd(   IA,JA,itp_nh       )
    integer :: kgrd(KA,IA,JA,itp_nh,itp_nv)

    real(RP) :: qc(KA,IA,JA)
    real(RP) :: qc_sfc(1,IA,JA)

    integer :: n, k, i, j, iq

    character(len=H_MID) :: varname_T
    character(len=H_MID) :: varname_W
    character(len=H_MID) :: varname_U
    character(len=H_MID) :: varname_V
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

    allocate( read_xy  (        dims(2),dims(3),ts:te   ) )
    allocate( read_uy  (        dims(5),dims(3),ts:te   ) )
    allocate( read_xv  (        dims(2),dims(6),ts:te   ) )
    allocate( read_zxy (dims(1),dims(2),dims(3),ts:te   ) )
    allocate( read_wxy (dims(4),dims(2),dims(3),ts:te   ) )
    allocate( read_zuy (dims(1),dims(5),dims(3),ts:te   ) )
    allocate( read_zxv (dims(1),dims(2),dims(6),ts:te   ) )

    allocate( velz_org (dims(4),dims(2),dims(3),ts:te   ) )
    allocate( velx_org (dims(1),dims(5),dims(3),ts:te   ) )
    allocate( vely_org (dims(1),dims(2),dims(6),ts:te   ) )
    allocate( pott_org (dims(1),dims(2),dims(3),ts:te   ) )
    allocate( temp_org (dims(1),dims(2),dims(3),ts:te   ) )
    allocate( pres_org (dims(1),dims(2),dims(3),ts:te   ) )
    allocate( p_org    (dims(1),dims(2),dims(3),ts:te   ) )
    allocate( pb_org   (dims(1),dims(2),dims(3),ts:te   ) )
    allocate( qtrc_org (dims(1),dims(2),dims(3),ts:te,QA) )

    allocate( tsfc_org (        dims(2),dims(3),ts:te   ) )
    allocate( psfc_org (        dims(2),dims(3),ts:te   ) )
    allocate( qsfc_org (        dims(2),dims(3),ts:te,QA) )

    allocate( ph_org   (dims(4),dims(2),dims(3),ts:te   ) )
    allocate( phb_org  (dims(4),dims(2),dims(3),ts:te   ) )
    allocate( geof_org (dims(4),dims(2),dims(3),ts:te   ) )
    allocate( geoh_org (dims(1),dims(2),dims(3),ts:te   ) )
    allocate( lat_org  (dims(2),dims(3),        ts:te   ) )
    allocate( lon_org  (dims(2),dims(3),        ts:te   ) )
    allocate( latu_org (dims(5),dims(3),        ts:te   ) )
    allocate( lonu_org (dims(5),dims(3),        ts:te   ) )
    allocate( latv_org (dims(2),dims(6),        ts:te   ) )
    allocate( lonv_org (dims(2),dims(6),        ts:te   ) )

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputWRF]'

    if( .NOT. serial .or. myrank == PRC_master ) then
       call ExternalFileRead( read_xy(:,:,:),    BASENAME, "XLAT",    ts, te, myrank, mdlid, single=.true.               )
       lat_org (:,:,:) = real( read_xy(:,:,:), kind=RP ) * D2R

       call ExternalFileRead( read_xy(:,:,:),    BASENAME, "XLONG",   ts, te, myrank, mdlid, single=.true.               )
       lon_org (:,:,:) = real( read_xy(:,:,:), kind=RP ) * D2R

       call ExternalFileRead( read_uy(:,:,:),    BASENAME, "XLAT_U",  ts, te, myrank, mdlid, single=.true., xstag=.true. )
       latu_org(:,:,:) = real( read_uy(:,:,:), kind=RP ) * D2R

       call ExternalFileRead( read_uy(:,:,:),    BASENAME, "XLONG_U", ts, te, myrank, mdlid, single=.true., xstag=.true. )
       lonu_org(:,:,:) = real( read_uy(:,:,:), kind=RP ) * D2R

       call ExternalFileRead( read_xv(:,:,:),    BASENAME, "XLAT_V",  ts, te, myrank, mdlid, single=.true., ystag=.true. )
       latv_org(:,:,:) = real( read_xv(:,:,:), kind=RP ) * D2R

       call ExternalFileRead( read_xv(:,:,:),    BASENAME, "XLONG_V", ts, te, myrank, mdlid, single=.true., ystag=.true. )
       lonv_org(:,:,:) = real( read_xv(:,:,:), kind=RP ) * D2R

       call ExternalFileRead( read_wxy(:,:,:,:), BASENAME, "PH",      ts, te, myrank, mdlid, single=.true., zstag=.true. )
       ph_org  (:,:,:,:) = real( read_wxy(:,:,:,:), kind=RP )

       call ExternalFileRead( read_wxy(:,:,:,:), BASENAME, "PHB",     ts, te, myrank, mdlid, single=.true., zstag=.true. )
       phb_org (:,:,:,:) = real( read_wxy(:,:,:,:), kind=RP )

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "P",       ts, te, myrank, mdlid, single=.true.               )
       p_org   (:,:,:,:) = real( read_zxy(:,:,:,:), kind=RP )

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "PB",      ts, te, myrank, mdlid, single=.true.               )
       pb_org  (:,:,:,:) = real( read_zxy(:,:,:,:), kind=RP )

       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, varname_T, ts, te, myrank, mdlid, single=.true.               )
       pott_org(:,:,:,:) = real( read_zxy(:,:,:,:), kind=RP )

       call ExternalFileRead( read_wxy(:,:,:,:), BASENAME, varname_W, ts, te, myrank, mdlid, single=.true., zstag=.true. )
       velz_org(:,:,:,:) = real( read_wxy(:,:,:,:), kind=RP )

       call ExternalFileRead( read_zuy(:,:,:,:), BASENAME, varname_U, ts, te, myrank, mdlid, single=.true., xstag=.true. )
       velx_org(:,:,:,:) = real( read_zuy(:,:,:,:), kind=RP )

       call ExternalFileRead( read_zxv(:,:,:,:), BASENAME, varname_V, ts, te, myrank, mdlid, single=.true., ystag=.true. )
       vely_org(:,:,:,:) = real( read_zxv(:,:,:,:), kind=RP )

       call ExternalFileRead( read_xy(:,:,:),    BASENAME, "PSFC",    ts, te, myrank, mdlid, single=.true.               )
       psfc_org (:,:,:) = real( read_xy(:,:,:), kind=RP )

       call ExternalFileRead( read_xy(:,:,:),    BASENAME, "T2",      ts, te, myrank, mdlid, single=.true.               )
       tsfc_org (:,:,:) = real( read_xy(:,:,:), kind=RP )

       call ExternalFileRead( read_xy(:,:,:),    BASENAME, "Q2",      ts, te, myrank, mdlid, single=.true.               )
       !convert from mixing ratio [kg/kg] to ratio of mass of tracer to total mass[kg/kg]
       qsfc_org (:,:,:,I_QV) = real(read_xy(:,:,:),kind=RP)/( 1.0_RP + real(read_xy(:,:,:),kind=RP) )

       qtrc_org(:,:,:,:,:) = 0.0_RP
       call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QVAPOR",  ts, te, myrank, mdlid, single=.true. )
       !convert from mixing ratio [kg/kg] to ratio of mass of tracer to total mass[kg/kg]
       qtrc_org(:,:,:,:,I_QV) = real(read_zxy(:,:,:,:),kind=RP)/( 1.0_RP + real(read_zxy(:,:,:,:),kind=RP) )

#ifndef DRY
       if( mptype_parent > 0 ) then
          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QCLOUD", ts, te, myrank, mdlid, single=.true. )
          !convert from mixing ratio [kg/kg] to ratio of mass of tracer to total mass[kg/kg]
          qtrc_org(:,:,:,:,I_QC) = real( read_zxy(:,:,:,:), kind=RP)/( 1.0_RP + real( read_zxy(:,:,:,:), kind=RP) )

          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QRAIN",  ts, te, myrank, mdlid, single=.true. )
          !convert from mixing ratio [kg/kg] to ratio of mass of tracer to total mass[kg/kg]
          qtrc_org(:,:,:,:,I_QR) = real( read_zxy(:,:,:,:), kind=RP)/( 1.0_RP + real( read_zxy(:,:,:,:), kind=RP) )
       endif

       if( mptype_parent > 3 ) then
          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QICE",   ts, te, myrank, mdlid, single=.true. )
          !convert from mixing ratio [kg/kg] to ratio of mass of tracer to total mass[kg/kg]
          qtrc_org(:,:,:,:,I_QI) = real( read_zxy(:,:,:,:), kind=RP)/( 1.0_RP + real( read_zxy(:,:,:,:), kind=RP) )

          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QSNOW",  ts, te, myrank, mdlid, single=.true. )
          !convert from mixing ratio [kg/kg] to ratio of mass of tracer to total mass[kg/kg]
          qtrc_org(:,:,:,:,I_QS) = real( read_zxy(:,:,:,:), kind=RP)/( 1.0_RP + real( read_zxy(:,:,:,:), kind=RP) )
       endif

       if( mptype_parent > 5 ) then
          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "QGRAUP", ts, te, myrank, mdlid, single=.true. )
          !convert from mixing ratio [kg/kg] to ratio of mass of tracer to total mass[kg/kg]
          qtrc_org(:,:,:,:,I_QG) = real( read_zxy(:,:,:,:), kind=RP)/( 1.0_RP + real( read_zxy(:,:,:,:), kind=RP) )
       endif

       if( mptype_parent > 6 ) then
          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NC",     ts, te, myrank, mdlid, single=.true. )
          qtrc_org(:,:,:,:,I_NC) = real( read_zxy(:,:,:,:), kind=RP )

          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NR",     ts, te, myrank, mdlid, single=.true. )
          qtrc_org(:,:,:,:,I_NR) = real( read_zxy(:,:,:,:), kind=RP )

          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NI",     ts, te, myrank, mdlid, single=.true. )
          qtrc_org(:,:,:,:,I_NI) = real( read_zxy(:,:,:,:), kind=RP )

          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NS",     ts, te, myrank, mdlid, single=.true. )
          qtrc_org(:,:,:,:,I_NS) = real( read_zxy(:,:,:,:), kind=RP )

          call ExternalFileRead( read_zxy(:,:,:,:), BASENAME, "NG",     ts, te, myrank, mdlid, single=.true. )
          qtrc_org(:,:,:,:,I_NG) = real( read_zxy(:,:,:,:), kind=RP )
       endif
#endif

       do n = ts, te
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)
          pres_org(k,i,j,n) = p_org(k,i,j,n) + pb_org(k,i,j,n)
          temp_org(k,i,j,n) = ( pott_org(k,i,j,n) + t0 ) * ( pres_org(k,i,j,n) / p0 )**RCP

          do iq = 1, QA
             qtrc_org(k,i,j,n,iq) = max( qtrc_org(k,i,j,n,iq), 0.0_RP )
          end do
       end do
       end do
       end do
       end do

       ! convert to geopotential height to use as real height in WRF
       do n = ts, te
       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(4)
          geof_org(k,i,j,n) = ( ph_org(k,i,j,n) + phb_org(k,i,j,n) ) / grav
       end do
       end do
       end do
       end do

       if( trim(mptype_run)=='double' .and. mptype_parent <= 6 )then
          if( IO_L ) write(IO_FID_LOG,*) '--- Diagnose Number Concentration from Mixing Ratio'
          call diagnose_number_concentration( qtrc_org(:,:,:,:,:) ) ! [inout]
       endif
    endif

    if( serial ) then
       call COMM_bcast( geof_org(:,:,:,:), dims(4),dims(2),dims(3),te )
       call COMM_bcast( lat_org (:,:,:),           dims(2),dims(3),te )
       call COMM_bcast( lon_org (:,:,:),           dims(2),dims(3),te )
       call COMM_bcast( latu_org(:,:,:),           dims(5),dims(3),te )
       call COMM_bcast( lonu_org(:,:,:),           dims(5),dims(3),te )
       call COMM_bcast( latv_org(:,:,:),           dims(2),dims(6),te )
       call COMM_bcast( lonv_org(:,:,:),           dims(2),dims(6),te )

       call COMM_bcast( velz_org(:,:,:,:), dims(4),dims(2),dims(3),te )
       call COMM_bcast( velx_org(:,:,:,:), dims(1),dims(5),dims(3),te )
       call COMM_bcast( vely_org(:,:,:,:), dims(1),dims(2),dims(6),te )
       call COMM_bcast( temp_org(:,:,:,:), dims(1),dims(2),dims(3),te )
       call COMM_bcast( pres_org(:,:,:,:), dims(1),dims(2),dims(3),te )
       call COMM_bcast( psfc_org(:,:,:),           dims(2),dims(3),te )
       call COMM_bcast( tsfc_org(:,:,:),           dims(2),dims(3),te )

       do iq = 1, QA
          call COMM_bcast( qtrc_org(:,:,:,:,iq), dims(1),dims(2),dims(3),te )
          call COMM_bcast( qsfc_org(:,:,:,iq),           dims(2),dims(3),te )
       enddo
    endif

    ! make half level of geopotential height from face level
    do n = ts, te
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       geoh_org(k,i,j,n) = ( geof_org(k,i,j,n) + geof_org(k+1,i,j,n) ) * 0.5_RP
    end do
    end do
    end do
    end do

    call check_domain_compatibility( lon_org(:,:,1),  lat_org(:,:,1),  geoh_org(:,:,:,1), &
                                     LON(:,:),        LAT(:,:),        CZ(:,:,:) )
    call check_domain_compatibility( lonu_org(:,:,1), latu_org(:,:,1), geoh_org(:,:,:,1), &
                                     LONX(:,:),       LAT(:,:),        CZ(:,:,:), skip_z=.true. )
    call check_domain_compatibility( lonv_org(:,:,1), latv_org(:,:,1), geoh_org(:,:,:,1), &
                                     LON(:,:),        LATY(:,:),       CZ(:,:,:), skip_z=.true. )
    call check_domain_compatibility( lon_org(:,:,1),  lat_org(:,:,1),  geof_org(:,:,:,1), &
                                     LON(:,:),        LAT(:,:),        FZ(:,:,:), skip_x=.true., skip_y=.true. )

    do n = ts, te !--- time loop

       ! for vector (w) points
       call latlonz_interpolation_fact( hfact   (:,:,:),           & ! [OUT]
                                        vfact   (:,:,:,:,:),       & ! [OUT]
                                        kgrd    (:,:,:,:,:),       & ! [OUT]
                                        igrd    (:,:,:),           & ! [OUT]
                                        jgrd    (:,:,:),           & ! [OUT]
                                        FZ      (:,:,:),           & ! [IN]
                                        LAT     (:,:),             & ! [IN]
                                        LON     (:,:),             & ! [IN]
                                        geof_org(:,:,:,:),         & ! [IN]
                                        lat_org (:,:,:),           & ! [IN]
                                        lon_org (:,:,:),           & ! [IN]
                                        dims(4), dims(2), dims(3), & ! [IN]
                                        n                          ) ! [IN]

       do j = 1, JA
       do i = 1, IA
       do k = KS, KE-1
          velz(k,i,j) = velz_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                      + velz_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                      + velz_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                      + velz_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                      + velz_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                      + velz_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
       end do
       end do
       end do

       ! for vector (u) points
       call latlonz_interpolation_fact( hfact   (:,:,:),           & ! [OUT]
                                        vfact   (:,:,:,:,:),       & ! [OUT]
                                        kgrd    (:,:,:,:,:),       & ! [OUT]
                                        igrd    (:,:,:),           & ! [OUT]
                                        jgrd    (:,:,:),           & ! [OUT]
                                        CZ      (:,:,:),           & ! [IN]
                                        LAT     (:,:),             & ! [IN]
                                        LON     (:,:),             & ! [IN]
                                        geoh_org(:,:,:,:),         & ! [IN]
                                        latu_org(:,:,:),           & ! [IN]
                                        lonu_org(:,:,:),           & ! [IN]
                                        dims(1), dims(2), dims(3), & ! [IN] ! dims(5) wasn't used to keep consistency with geoh_org
                                        n                          ) ! [IN]

       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          velx(k,i,j) = velx_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                      + velx_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                      + velx_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                      + velx_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                      + velx_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                      + velx_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
       end do
       end do
       end do

       ! for vector (v) points
       call latlonz_interpolation_fact( hfact   (:,:,:),           & ! [OUT]
                                        vfact   (:,:,:,:,:),       & ! [OUT]
                                        kgrd    (:,:,:,:,:),       & ! [OUT]
                                        igrd    (:,:,:),           & ! [OUT]
                                        jgrd    (:,:,:),           & ! [OUT]
                                        CZ      (:,:,:),           & ! [IN]
                                        LAT     (:,:),             & ! [IN]
                                        LON     (:,:),             & ! [IN]
                                        geoh_org(:,:,:,:),         & ! [IN]
                                        latv_org(:,:,:),           & ! [IN]
                                        lonv_org(:,:,:),           & ! [IN]
                                        dims(1), dims(2), dims(3), & ! [IN] ! dims(6) wasn't used to keep consistency with geoh_org
                                        n                          ) ! [IN]

       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          vely(k,i,j) = vely_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                      + vely_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                      + vely_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                      + vely_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                      + vely_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                      + vely_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
       end do
       end do
       end do

       call wrf_arwpost_calc_uvmet( llvelx, llvely, & ! (out)
                                    velx,   vely,   & ! (in)
                                    LON, LAT,       & ! (in)
                                    BASENAME        ) ! (in)

       ! convert from latlon coordinate to local mapping (x)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          work(k,i,j) = llvelx(k,i,j) * rotc(i,j,cosin) + llvely(k,i,j) * rotc(i,j,sine)
       end do
       end do
       end do

       ! from scalar point to staggered point
       do j = 1, JA
       do i = 1, IA-1
       do k = KS, KE
          velx(k,i,j) = ( work(k,i+1,j) + work(k,i,j) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do k = KS, KE
          velx(k,IA,j) = work(k,IA,j)
       end do
       end do
       velx(KS-1,:,:) = 0.0_RP
       velx(KS-2,:,:) = 0.0_RP
       call COMM_vars8( velx(:,:,:), 1 )

       ! convert from latlon coordinate to local mapping (y)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          work(k,i,j) = - llvelx(k,i,j) * rotc(i,j,sine) + llvely(k,i,j) * rotc(i,j,cosin)
       end do
       end do
       end do

       ! from scalar point to staggered point
       do j = 1, JA-1
       do i = 1, IA
       do k = KS, KE
          vely(k,i,j) = ( work(k,i,j+1) + work(k,i,j) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          vely(k,i,JA) = work(k,i,JA)
       end do
       end do
       vely(KS-1,:,:) = 0.0_RP
       vely(KS-2,:,:) = 0.0_RP
       call COMM_vars8( vely(:,:,:), 2 )

       call COMM_wait ( velx(:,:,:), 1, .false. )
       call COMM_wait ( vely(:,:,:), 2, .false. )

       ! for scalar points
       call latlonz_interpolation_fact( hfact   (:,:,:),           & ! [OUT]
                                        vfact   (:,:,:,:,:),       & ! [OUT]
                                        kgrd    (:,:,:,:,:),       & ! [OUT]
                                        igrd    (:,:,:),           & ! [OUT]
                                        jgrd    (:,:,:),           & ! [OUT]
                                        CZ      (:,:,:),           & ! [IN]
                                        LAT     (:,:),             & ! [IN]
                                        LON     (:,:),             & ! [IN]
                                        geoh_org(:,:,:,:),         & ! [IN]
                                        lat_org (:,:,:),           & ! [IN]
                                        lon_org (:,:,:),           & ! [IN]
                                        dims(1), dims(2), dims(3), & ! [IN]
                                        n                          ) ! [IN]

       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          temp(k,i,j) = temp_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                      + temp_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                      + temp_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                      + temp_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                      + temp_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                      + temp_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)

          pres(k,i,j) = pres_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                      + pres_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                      + pres_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                      + pres_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                      + pres_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                      + pres_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)

          do iq = 1, QA
             qtrc(k,i,j,n,iq) = qtrc_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n,iq) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                              + qtrc_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n,iq) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                              + qtrc_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n,iq) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                              + qtrc_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n,iq) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                              + qtrc_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n,iq) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                              + qtrc_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n,iq) * hfact(i,j,3) * vfact(k,i,j,3,2)
          end do

          call THERMODYN_pott( pott(k,i,j),    & ! [OUT]
                               temp(k,i,j),    & ! [IN]
                               pres(k,i,j),    & ! [IN]
                               qtrc(k,i,j,n,:) ) ! [IN]
       end do
       end do
       end do

       do j = 1, JA
       do i = 1, IA
          temp_sfc(1,i,j) = tsfc_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                          + tsfc_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                          + tsfc_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)

          pres_sfc(1,i,j) = psfc_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                          + psfc_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                          + psfc_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)

          ! interpolate QV (=Q2) only: other QTRC are set zero
          qtrc_sfc(1,i,j,:)    = 0.0_RP
          qtrc_sfc(1,i,j,I_QV) = qsfc_org(igrd(i,j,1),jgrd(i,j,1),n,I_QV) * hfact(i,j,1) &
                               + qsfc_org(igrd(i,j,2),jgrd(i,j,2),n,I_QV) * hfact(i,j,2) &
                               + qsfc_org(igrd(i,j,3),jgrd(i,j,3),n,I_QV) * hfact(i,j,3)

          call THERMODYN_pott( pott_sfc(1,i,j),  & ! [OUT]
                               temp_sfc(1,i,j),  & ! [IN]
                               pres_sfc(1,i,j),  & ! [IN]
                               qtrc_sfc(1,i,j,:) ) ! [IN]
       end do
       end do

#ifdef DRY
       qc = 0.0_RP
       qc_sfc = 0.0_RP
#else
       qc = qtrc(:,:,:,n,I_QC)
       qc_sfc = qtrc_sfc(:,:,:,I_QC)
#endif
       ! make density & pressure profile in moist condition
       call HYDROSTATIC_buildrho_real( dens    (:,:,:,n),      & ! [OUT]
                                       temp    (:,:,:),        & ! [OUT]
                                       pres    (:,:,:),        & ! [OUT]
                                       pott    (:,:,:),        & ! [IN]
                                       qtrc    (:,:,:,n,I_QV), & ! [IN]
                                       qc      (:,:,:),        & ! [IN]
                                       temp_sfc(:,:,:),        & ! [OUT]
                                       pres_sfc(:,:,:),        & ! [IN]
                                       pott_sfc(:,:,:),        & ! [IN]
                                       qtrc_sfc(:,:,:,I_QV),   & ! [IN]
                                       qc_sfc  (:,:,:)         ) ! [IN]

       call COMM_vars8( dens(:,:,:,n), 1 )
       call COMM_wait ( dens(:,:,:,n), 1 )

       do j = 1, JA
       do i = 1, IA
       do k = KS, KE-1
          momz(k,i,j,n) = 0.5_RP * velz(k,i,j) * ( dens(k+1,  i,  j,n) + dens(k,i,j,n) )
       end do
       end do
       end do
       do j = 1, JA
       do i = 1, IA
          momz(KE,i,j,n) = 0.0_RP
       end do
       end do
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          rhot(k,i,j,n) = pott(k,i,j) * dens(k,i,j,n)
       end do
       end do
       end do

       do j = 1, JA
       do i = 1, IA-1
       do k = KS, KE
          momx(k,i,j,n) = 0.5_RP * velx(k,i,j) * ( dens(k,  i+1,  j,n) + dens(k,i,j,n) )
       end do
       end do
       end do
       do j = 1, JA
       do k = KS, KE
          momx(k,IA,j,n) = velx(k,IA,j) * dens(k,IA,j,n)
       end do
       end do
       call COMM_vars8( momx(:,:,:,n), 1 )

       do j = 1, JA-1
       do i = 1, IA
       do k = KS, KE
          momy(k,i,j,n) = 0.5_RP * vely(k,i,j) * ( dens(k,    i,j+1,n) + dens(k,i,j,n) )
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          momy(k,i,JA,n) = vely(k,i,JA) * dens(k,i,JA,n)
       end do
       end do
       call COMM_vars8( momy(:,:,:,n), 2 )

       call COMM_wait ( momx(:,:,:,n), 1, .false. )
       call COMM_wait ( momy(:,:,:,n), 2, .false. )

    end do !--- time loop

    deallocate( read_xy  )
    deallocate( read_uy  )
    deallocate( read_xv  )
    deallocate( read_zxy )
    deallocate( read_wxy )
    deallocate( read_zuy )
    deallocate( read_zxv )

    deallocate( velz_org )
    deallocate( velx_org )
    deallocate( vely_org )
    deallocate( pott_org )
    deallocate( temp_org )
    deallocate( pres_org )
    deallocate( p_org    )
    deallocate( pb_org   )
    deallocate( qtrc_org )

    deallocate( ph_org   )
    deallocate( phb_org  )
    deallocate( geof_org )
    deallocate( geoh_org )
    deallocate( lat_org  )
    deallocate( lon_org  )
    deallocate( latu_org )
    deallocate( lonu_org )
    deallocate( latv_org )
    deallocate( lonv_org )

    return
  end subroutine InputAtomWRF


  !-----------------------------------------------------------------------------
  !> Individual Procedures
  !-----------------------------------------------------------------------------
  subroutine InputAtomNICAM( &
      dens,             & ! (out)
      momz,             & ! (out)
      momx,             & ! (out)
      momy,             & ! (out)
      rhot,             & ! (out)
      qtrc,             & ! (out)
      basename_num,     & ! (in)
      dims,             & ! (in)
      start_step,       & ! (in)
      end_step          ) ! (in)
    use scale_const, only: &
       D2R => CONST_D2R,   &
       EPS => CONST_EPS
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
       THERMODYN_pott      => ATMOS_THERMODYN_pott
    use scale_gridtrans, only: &
       rotc => GTRANS_ROTC
    use scale_topography, only: &
       topo => TOPO_Zsfc
    implicit none

    real(RP),         intent(out) :: dens(:,:,:,:)
    real(RP),         intent(out) :: momz(:,:,:,:)
    real(RP),         intent(out) :: momx(:,:,:,:)
    real(RP),         intent(out) :: momy(:,:,:,:)
    real(RP),         intent(out) :: rhot(:,:,:,:)
    real(RP),         intent(out) :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in)  :: basename_num
    integer,          intent(in)  :: dims(:)
    integer,          intent(in)  :: start_step
    integer,          intent(in)  :: end_step

    !> work
    real(SP), allocatable :: read1DX(:)
    real(SP), allocatable :: read1DY(:)
    real(SP), allocatable :: read1DZ(:)
    real(SP), allocatable :: read3DT(:,:,:,:)
    real(SP), allocatable :: read4D (:,:,:,:)

    real(RP), allocatable :: lon_org (:,:,:)
    real(RP), allocatable :: lat_org (:,:,:)
    real(RP), allocatable :: hgt_org(:,:,:,:)
    real(RP), allocatable :: hgt_op1(:,:,:,:)

    real(RP), allocatable :: tsfc_org(:,:,:)
    real(RP), allocatable :: slp_org(:,:,:)
    !real(RP), allocatable :: psfc_org(:,:,:)
    real(RP), allocatable :: qsfc_org(:,:,:,:)

    real(RP), allocatable :: velx_org(:,:,:,:)
    real(RP), allocatable :: vely_org(:,:,:,:)
    real(RP), allocatable :: temp_org(:,:,:,:)
    real(RP), allocatable :: pres_org(:,:,:,:)
    real(RP), allocatable :: qtrc_org(:,:,:,:,:)

    real(RP) :: velz  (KA,IA,JA,start_step:end_step)
    real(RP) :: velx  (KA,IA,JA,start_step:end_step)
    real(RP) :: vely  (KA,IA,JA,start_step:end_step)
    real(RP) :: llvelx(KA,IA,JA,start_step:end_step)
    real(RP) :: llvely(KA,IA,JA,start_step:end_step)
    real(RP) :: pott  (KA,IA,JA,start_step:end_step)
    real(RP) :: temp  (KA,IA,JA,start_step:end_step)
    real(RP) :: pres  (KA,IA,JA,start_step:end_step)
    real(RP) :: work  (KA,IA,JA,start_step:end_step)

    real(RP) :: pott_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: slp_sfc (1,IA,JA,start_step:end_step)
    real(RP) :: pres_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: temp_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: qtrc_sfc(1,IA,JA,start_step:end_step,QA)

    real(RP) :: hfact(   IA,JA,itp_nh       )
    real(RP) :: vfact(KA,IA,JA,itp_nh,itp_nv)

    integer :: igrd(   IA,JA,itp_nh       )
    integer :: jgrd(   IA,JA,itp_nh       )
    integer :: kgrd(KA,IA,JA,itp_nh,itp_nv)

    real(RP) :: qc(KA,IA,JA)
    real(RP) :: qc_sfc(1,IA,JA)

    real(RP) :: sw
    real(RP) :: Rtot
    real(RP) :: wgt_up, wgt_bt
    real(RP) :: z1,z2,pres1,pres2

    integer :: step_fixed = 1
    integer :: bottom, upper
    integer :: QA_NCM = 1
    integer :: dims3p1

    integer :: k, i, j, n
    integer :: kk, kks, kke
    integer :: iq, ierr
    integer :: nt
    logical :: do_read, lack_of_val

    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    nt = end_step - start_step + 1
    dims3p1 = dims(3) + 1

    if( myrank == PRC_master ) then
       do_read = .true.
    else
       do_read = .false.
    endif

    if( do_read ) then
       allocate( read1DX( dims(1)                       ) )
       allocate( read1DY( dims(2)                       ) )
       allocate( read1DZ( dims(3)                       ) )
       allocate( read3DT( 1,       dims(1), dims(2), nt ) )
       allocate( read4D ( dims(3), dims(1), dims(2), nt ) )
    endif

    allocate( lon_org(           dims(1), dims(2), step_fixed ) )
    allocate( lat_org(           dims(1), dims(2), step_fixed ) )
    allocate( hgt_org( dims(3),  dims(1), dims(2), step_fixed ) )
    allocate( hgt_op1( dims3p1,  dims(1), dims(2), step_fixed ) )

    allocate( tsfc_org(          dims(1), dims(2), start_step:end_step ) )
    allocate( slp_org (          dims(1), dims(2), start_step:end_step ) )
    !allocate( psfc_org(          dims(1), dims(2), start_step:end_step ) )
    allocate( qsfc_org(          dims(1), dims(2), start_step:end_step, QA_NCM) )
    allocate( velx_org( dims(3), dims(1), dims(2), start_step:end_step ) )
    allocate( vely_org( dims(3), dims(1), dims(2), start_step:end_step ) )
    allocate( temp_org( dims3p1, dims(1), dims(2), start_step:end_step ) )
    allocate( pres_org( dims3p1, dims(1), dims(2), start_step:end_step ) )
    allocate( qtrc_org( dims(3), dims(1), dims(2), start_step:end_step, QA_NCM) )

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputNICAM]'

    ! normal vertical grid arrangement
    if( do_read ) then
       basename = "ms_pres"//trim(basename_num)
       call FileRead( read1DX(:), trim(basename), "lon", step_fixed, 1, single=.true. )
       do j = 1, dims(2)
          lon_org (:,j,step_fixed)  = read1DX(:) * D2R
       enddo

       call FileRead( read1DY(:), trim(basename), "lat", step_fixed, 1, single=.true. )
       do i = 1, dims(1)
          lat_org (i,:,step_fixed)  = read1DY(:) * D2R
       enddo

       call FileRead( read1DZ(:), trim(basename), "lev", step_fixed, 1, single=.true. )
       do j = 1, dims(2)
       do i = 1, dims(1)
          hgt_org(:,i,j,step_fixed) = read1DZ(:)
       end do
       end do
    endif
    call COMM_bcast( lat_org (:,:,:),                 dims(1), dims(2), step_fixed )
    call COMM_bcast( lon_org (:,:,:),                 dims(1), dims(2), step_fixed )
    call COMM_bcast( hgt_org (:,:,:,:),      dims(3), dims(1), dims(2), step_fixed )

    call check_domain_compatibility( lon_org(:,:,1),  lat_org(:,:,1),  hgt_org(:,:,:,1), &
                                     LON(:,:),        LAT(:,:),        CZ(:,:,:) )
    call check_domain_compatibility( lon_org(:,:,1),  lat_org(:,:,1),  hgt_org(:,:,:,1), &
                                     LONX(:,:),       LAT(:,:),        CZ(:,:,:), skip_y=.true., skip_z=.true. )
    call check_domain_compatibility( lon_org(:,:,1),  lat_org(:,:,1),  hgt_org(:,:,:,1), &
                                     LON(:,:),        LATY(:,:),       CZ(:,:,:), skip_x=.true., skip_z=.true. )
    call check_domain_compatibility( lon_org(:,:,1),  lat_org(:,:,1),  hgt_org(:,:,:,1), &
                                     LON(:,:),        LAT(:,:),        FZ(:,:,:), skip_x=.true., skip_y=.true. )

    call latlonz_interpolation_fact( hfact   (:,:,:),           & ! [OUT]
                                     vfact   (:,:,:,:,:),       & ! [OUT]
                                     kgrd    (:,:,:,:,:),       & ! [OUT]
                                     igrd    (:,:,:),           & ! [OUT]
                                     jgrd    (:,:,:),           & ! [OUT]
                                     CZ      (:,:,:),           & ! [IN]
                                     LAT     (:,:),             & ! [IN]
                                     LON     (:,:),             & ! [IN]
                                     hgt_org (:,:,:,:),         & ! [IN]
                                     lat_org (:,:,:),           & ! [IN]
                                     lon_org (:,:,:),           & ! [IN]
                                     dims(3), dims(1), dims(2), & ! [IN]
                                     step_fixed                 ) ! [IN]
    deallocate( hgt_org  )

    if( do_read ) then
       !> [scale-offset]
       basename = "ms_u"//trim(basename_num)
       call ExternalFileReadOffset( read4D(:,:,:,:),  &
                                    trim(basename),   &
                                    "ms_u",           &
                                    start_step,       &
                                    end_step,         &
                                    myrank,           &
                                    iNICAM,           &
                                    single=.true.     )
       velx_org(:,:,:,:) = real( read4D(:,:,:,:), kind=RP )
    endif
    call COMM_bcast( velx_org(:,:,:,:),      dims(3), dims(1), dims(2), nt )

    ! tentative: missing value => 0.
    do n = start_step, end_step
    do j = 1, dims(2)
    do i = 1, dims(1)
    do k = 1, dims(3)
       if( abs(abs(velx_org( k, i, j, n ))-300.0D0) < sqrt(EPS) )then
           velx_org( k, i, j, n ) = 0.0D0
       endif
    end do
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
       llvelx(k,i,j,n) = velx_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                       + velx_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                       + velx_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                       + velx_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                       + velx_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                       + velx_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    end do
    end do
    end do
    end do
    deallocate( velx_org )

    if( do_read ) then
       !> [scale-offset]
       basename = "ms_v"//trim(basename_num)
       call ExternalFileReadOffset( read4D(:,:,:,:),  &
                                    trim(basename),   &
                                    "ms_v",           &
                                    start_step,       &
                                    end_step,         &
                                    myrank,           &
                                    iNICAM,           &
                                    single=.true.     )
       vely_org(:,:,:,:) = real( read4D(:,:,:,:), kind=RP )
    endif
    call COMM_bcast( vely_org(:,:,:,:),      dims(3), dims(1), dims(2), nt )

    ! tentative: missing value => 0.
    do n = start_step, end_step
    do j = 1, dims(2)
    do i = 1, dims(1)
    do k = 1, dims(3)
       if( abs(abs(vely_org( k, i, j, n ))-300.0D0) < sqrt(EPS) )then
           vely_org( k, i, j, n ) = 0.0D0
       endif
    end do
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
       llvely(k,i,j,n) = vely_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                       + vely_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                       + vely_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                       + vely_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                       + vely_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                       + vely_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    end do
    end do
    end do
    end do
    deallocate( vely_org )

    do n = start_step, end_step
       ! convert from latlon coordinate to local mapping (x)
       do j = 1, JA
       do i = 1, IA
       do k = KS-1, KE+1
          work(k,i,j,n) = llvelx(k,i,j,n) * rotc(i,j,cosin) + llvely(k,i,j,n) * rotc(i,j,sine)
       end do
       end do
       end do

       ! from scalar point to staggered point
       do j = 1, JA
       do i = 1, IA-1
       do k = KS, KE
          velx(k,i,j,n) = ( work(k,i+1,j,n) + work(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do k = KS, KE
          velx(k,IA,j,n) = work(k,IA,j,n)
       end do
       end do

       velx(KS-1,:,:,n) = 0.0_RP
       velx(KS-2,:,:,n) = 0.0_RP

       call COMM_vars8( velx(:,:,:,n), 1 )
       call COMM_wait ( velx(:,:,:,n), 1, .false. )

    end do

    do n = start_step, end_step
       ! convert from latlon coordinate to local mapping (y)
       do j = 1, JA
       do i = 1, IA
       do k = KS-1, KE+1
          work(k,i,j,n) = - llvelx(k,i,j,n) * rotc(i,j,sine) + llvely(k,i,j,n) * rotc(i,j,cosin)
       end do
       end do
       end do

       ! from scalar point to staggered point
       do j = 1, JA-1
       do i = 1, IA
       do k = KS, KE
          vely(k,i,j,n) = ( work(k,i,j+1,n) + work(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          vely(k,i,JA,n) = work(k,i,JA,n)
       end do
       end do
       vely(KS-1,:,:,n) = 0.0_RP
       vely(KS-2,:,:,n) = 0.0_RP

       call COMM_vars8( vely(:,:,:,n), 1 )
       call COMM_wait ( vely(:,:,:,n), 1, .false. )

    end do

    velz(:,:,:,:)   = 0.0_RP !> cold initialize for vertical velocity

    if( do_read ) then
       !> [scale-offset]
       basename = "ss_slp"//trim(basename_num)
       call ExternalFileReadOffset( read3DT(:,:,:,:), &
                                    trim(basename),   &
                                    "ss_slp",         &
                                    start_step,       &
                                    end_step,         &
                                    myrank,           &
                                    iNICAM,           &
                                    single=.true.     )
       slp_org(:,:,:) = real( read3DT(1,:,:,:), kind=RP )

       basename = "ss_t2m"//trim(basename_num)
       call ExternalFileRead( read3DT(:,:,:,:), trim(basename), "ss_t2m", start_step, end_step, myrank, iNICAM, single=.true. )
       tsfc_org(:,:,:) = real( read3DT(1,:,:,:), kind=RP )

       basename = "ss_q2m"//trim(basename_num)
       call ExternalFileRead( read3DT(:,:,:,:), trim(basename), "ss_q2m", start_step, end_step, myrank, iNICAM, single=.true. )
       qsfc_org(:,:,:,I_QV) = real( read3DT(1,:,:,:), kind=RP )
    endif
    call COMM_bcast( slp_org(:,:,:),                  dims(1), dims(2), nt )
    call COMM_bcast( tsfc_org(:,:,:),                 dims(1), dims(2), nt )
    call COMM_bcast( qsfc_org(:,:,:,I_QV),            dims(1), dims(2), nt )

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
       !pres_sfc(1,i,j,n) = psfc_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
       !                  + psfc_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
       !                  + psfc_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)

       slp_sfc(1,i,j,n) =  slp_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                         + slp_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                         + slp_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)

       temp_sfc(1,i,j,n) = tsfc_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                         + tsfc_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                         + tsfc_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)

       !> QTRC_sfc are set zero, except for QV; set q2.
       qtrc_sfc(1,i,j,n,:)    = 0.0_RP
       qtrc_sfc(1,i,j,n,I_QV) = qsfc_org(igrd(i,j,1),jgrd(i,j,1),n,I_QV) * hfact(i,j,1) &
                              + qsfc_org(igrd(i,j,2),jgrd(i,j,2),n,I_QV) * hfact(i,j,2) &
                              + qsfc_org(igrd(i,j,3),jgrd(i,j,3),n,I_QV) * hfact(i,j,3)
       sw = sign(0.5_RP, qtrc_sfc(1,i,j,n,I_QV)) + 0.5_RP
       qtrc_sfc(1,i,j,n,I_QV) = qtrc_sfc(1,i,j,n,I_QV) * sw !> --fix negative value

       !call THERMODYN_pott( pott_sfc(1,i,j,n),  & ! [OUT]
       !                     temp_sfc(1,i,j,n),  & ! [IN]
       !                     pres_sfc(1,i,j,n),  & ! [IN]
       !                     qtrc_sfc(1,i,j,n,:) ) ! [IN]
    end do
    end do
    end do


    if( do_read ) then
       basename = "ms_qv"//trim(basename_num)
       call ExternalFileRead( read4D(:,:,:,:), trim(basename), "ms_qv", start_step, end_step, myrank, iNICAM, single=.true. )
       qtrc_org(:,:,:,:,I_QV) = real( read4D(:,:,:,:), kind=RP )
    endif
    call COMM_bcast( qtrc_org(:,:,:,:,I_QV), dims(3), dims(1), dims(2), nt )

       ! tentative: missing value => surface_qv
       do n = start_step, end_step
       do j = 1, dims(2)
       do i = 1, dims(1)
       do k = 1, dims(3)
          if( qtrc_org(k,i,j,n,I_QV) < -9.99*10.**10 )then
             qtrc_org(k,i,j,n,I_QV) = qsfc_org(i,j,n,I_QV)
          endif
       end do
       end do
       end do
       end do

       do n = start_step, end_step
       do j = 1, dims(2)
       do i = 1, dims(1)
       do k = 1, dims(3)
          !> --fix negative value
          sw = sign(0.5_RP, qtrc_org(k,i,j,n,I_QV)) + 0.5_RP
          qtrc_org(k,i,j,n,I_QV) = qtrc_org(k,i,j,n,I_QV) * sw
       enddo
       enddo
       enddo
       enddo

    qtrc(:,:,:,:,:) = 0.0_RP !> cold initialize for hydrometeor
    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
       qtrc(k,i,j,n,I_QV) = qtrc_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n,I_QV) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                          + qtrc_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n,I_QV) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                          + qtrc_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n,I_QV) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                          + qtrc_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n,I_QV) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                          + qtrc_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n,I_QV) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                          + qtrc_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n,I_QV) * hfact(i,j,3) * vfact(k,i,j,3,2)
    end do
    end do
    end do
    end do
    deallocate( qtrc_org )

    ! vertical grid arrangement combined with surface
    if( do_read ) then
       do j = 1, dims(2)
       do i = 1, dims(1)
          hgt_op1(1,i,j,step_fixed) = 0.0_RP !bottom level is 0 [m]
          do k = 2, dims3p1
             hgt_op1(k,i,j,step_fixed) = read1DZ(k-1)
          enddo
       end do
       end do
    endif
    call COMM_bcast( hgt_op1(:,:,:,:), dims3p1, dims(1), dims(2), step_fixed )

    call latlonz_interpolation_fact( hfact   (:,:,:),           & ! [OUT]
                                     vfact   (:,:,:,:,:),       & ! [OUT]
                                     kgrd    (:,:,:,:,:),       & ! [OUT]
                                     igrd    (:,:,:),           & ! [OUT]
                                     jgrd    (:,:,:),           & ! [OUT]
                                     CZ      (:,:,:),           & ! [IN]
                                     LAT     (:,:),             & ! [IN]
                                     LON     (:,:),             & ! [IN]
                                     hgt_op1 (:,:,:,:),         & ! [IN]
                                     lat_org (:,:,:),           & ! [IN]
                                     lon_org (:,:,:),           & ! [IN]
                                     dims3p1, dims(1), dims(2), & ! [IN]
                                     step_fixed                 ) ! [IN]

    if( do_read ) then
       basename = "ms_pres"//trim(basename_num)
       call ExternalFileRead( read4D(:,:,:,:), trim(basename), "ms_pres", start_step, end_step, myrank, iNICAM, single=.true. )
       pres_org(1,:,:,:) = slp_org(:,:,:)
       do n = start_step, end_step
       do j = 1, dims(2)
       do i = 1, dims(1)
       do k = 2, dims3p1
          pres_org(k,i,j,n) = real( read4D(k-1,i,j,n), kind=RP )
       enddo
       enddo
       enddo
       enddo

       ! interpolate using SLP (fill lacked data)
       do n = start_step, end_step
       do j = 1, dims(2)
       do i = 1, dims(1)
          lack_of_val = .false.

          do k = 2, dims3p1
             if ( pres_org(k,i,j,n) <= 0.0_RP ) then
                if ( .NOT. lack_of_val ) kks = k
                kke = k
                lack_of_val = .true.
             endif
          enddo

          if ( lack_of_val ) then
             bottom = 1
             upper  = kke + 1
             do k = kks, kke
                wgt_bt = (hgt_op1(upper,i,j,1) - hgt_op1(k,     i,j,1)) / ((hgt_op1(upper,i,j,1) - hgt_op1(bottom,i,j,1)))
                wgt_up = (hgt_op1(k,    i,j,1) - hgt_op1(bottom,i,j,1)) / ((hgt_op1(upper,i,j,1) - hgt_op1(bottom,i,j,1)))
                pres_org(k,i,j,n) = exp( log(pres_org(bottom,i,j,n))*wgt_bt + log(pres_org(kke+1,i,j,n))*wgt_up )
             enddo
          endif
       enddo
       enddo
       enddo
    endif
    call COMM_bcast( pres_org(:,:,:,:), dims3p1, dims(1), dims(2), nt )

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
       pres(k,i,j,n) = pres_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + pres_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + pres_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + pres_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + pres_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + pres_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    end do
    end do
    end do
    end do
    deallocate( pres_org )

    ! Interpolate Surface pressure from SLP and PRES
    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
       lack_of_val = .true.
       do k = KS-1, KE
          if(k == KS-1)then
            z1=0.0_RP
            z2=CZ(k+1,i,j)
            pres1=slp_sfc(1,i,j,n)
            pres2=pres(k+1,i,j,n)
          else
            z1=CZ(k,i,j)
            z2=CZ(k+1,i,j)
            pres1=pres(k,i,j,n)
            pres2=pres(k+1,i,j,n)
          endif
          if((topo(i,j)>=z1).and.(topo(i,j)<z2))then
             lack_of_val = .false.                  ! found
             wgt_bt = (z2        - topo(i,j)) / (z2 - z1)
             wgt_up = (topo(i,j) - z1       ) / (z2 - z1)
             pres_sfc(1,i,j,n) = exp( log(pres1)*wgt_bt + log(pres2)*wgt_up )
             exit ! exit loop
          endif
       enddo
       if( lack_of_val )then
          write(IO_FID_LOG,*) 'realinput atom NICAM : cannot estimate pres_sfc',i,j,n
          call PRC_MPIstop
       endif
    enddo
    enddo
    enddo

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
       call THERMODYN_pott( pott_sfc(1,i,j,n),  & ! [OUT]
                            temp_sfc(1,i,j,n),  & ! [IN]
                            pres_sfc(1,i,j,n),  & ! [IN]
                            qtrc_sfc(1,i,j,n,:) ) ! [IN]
    end do
    end do
    end do


    ! vertical grid arrangement combined with 2m height
    if( do_read ) then
       do j = 1, dims(2)
       do i = 1, dims(1)
          hgt_op1(1,i,j,step_fixed) = 2.0_RP !bottom level is 2 [m]
          do k = 2, dims3p1
             hgt_op1(k,i,j,step_fixed) = read1DZ(k-1)
          enddo
       end do
       end do
    endif
    call COMM_bcast( hgt_op1(:,:,:,:), dims3p1, dims(1), dims(2), step_fixed )

    call latlonz_interpolation_fact( hfact   (:,:,:),           & ! [OUT]
                                     vfact   (:,:,:,:,:),       & ! [OUT]
                                     kgrd    (:,:,:,:,:),       & ! [OUT]
                                     igrd    (:,:,:),           & ! [OUT]
                                     jgrd    (:,:,:),           & ! [OUT]
                                     CZ      (:,:,:),           & ! [IN]
                                     LAT     (:,:),             & ! [IN]
                                     LON     (:,:),             & ! [IN]
                                     hgt_op1 (:,:,:,:),         & ! [IN]
                                     lat_org (:,:,:),           & ! [IN]
                                     lon_org (:,:,:),           & ! [IN]
                                     dims3p1, dims(1), dims(2), & ! [IN]
                                     step_fixed                 ) ! [IN]
    deallocate( lon_org  )
    deallocate( lat_org  )

    if( do_read ) then
       !> [scale-offset]
       basename = "ms_tem"//trim(basename_num)
       call ExternalFileReadOffset( read4D(:,:,:,:),  &
                                    trim(basename),   &
                                    "ms_tem",         &
                                    start_step,       &
                                    end_step,         &
                                    myrank,           &
                                    iNICAM,           &
                                    single=.true.     )
       temp_org(1,:,:,:) = tsfc_org(:,:,:)
       do n = start_step, end_step
       do j = 1, dims(2)
       do i = 1, dims(1)
       do k = 2, dims3p1
          temp_org(k,i,j,n) = real( read4D(k-1,i,j,n), kind=RP )
       enddo
       enddo
       enddo
       enddo

       ! interpolate using T2m (fill lacked data)
       do n = start_step, end_step
       do j = 1, dims(2)
       do i = 1, dims(1)
          lack_of_val = .false.

          do k = 2, dims3p1
             if ( temp_org(k,i,j,n) <= 50.0_RP ) then
                if ( .NOT. lack_of_val ) kks = k
                kke = k
                lack_of_val = .true.
             endif
          enddo

          if ( lack_of_val ) then
             bottom = 1
             upper  = kke + 1
             do k = kks, kke
                wgt_bt = (hgt_op1(upper,i,j,1) - hgt_op1(k,     i,j,1)) / ((hgt_op1(upper,i,j,1) - hgt_op1(bottom,i,j,1)))
                wgt_up = (hgt_op1(k,    i,j,1) - hgt_op1(bottom,i,j,1)) / ((hgt_op1(upper,i,j,1) - hgt_op1(bottom,i,j,1)))
                temp_org(k,i,j,n) = temp_org(bottom,i,j,n)*wgt_bt + temp_org(kke+1,i,j,n)*wgt_up
             enddo
          endif
       enddo
       enddo
       enddo
    endif
    call COMM_bcast( temp_org(:,:,:,:), dims3p1, dims(1), dims(2), nt )

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
       temp(k,i,j,n) = temp_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + temp_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + temp_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + temp_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + temp_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + temp_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    end do
    end do
    end do
    end do
    deallocate( temp_org )
    deallocate( slp_org )
    deallocate( tsfc_org )
    deallocate( qsfc_org )
    deallocate( hgt_op1  )

    if( do_read ) then
       deallocate( read1DX )
       deallocate( read1DY )
       deallocate( read1DZ )
       deallocate( read3DT )
       deallocate( read4D  )
    endif

    do n = start_step, end_step
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE+1
          call THERMODYN_pott( pott(k,i,j,n),  & ! [OUT]
                               temp(k,i,j,n),  & ! [IN]
                               pres(k,i,j,n),  & ! [IN]
                               qtrc(k,i,j,n,:) ) ! [IN]
    end do
    end do
    end do
    end do

    do n = start_step, end_step
#ifdef DRY
       qc = 0.0_RP
       qc_sfc = 0.0_RP
#else
       qc = qtrc(:,:,:,n,I_QC)
       qc_sfc = qtrc_sfc(:,:,:,n,I_QC)
#endif
       !> make density in moist condition
       call HYDROSTATIC_buildrho_real( dens    (:,:,:,n),      & ! [OUT]
                                       temp    (:,:,:,n),      & ! [OUT] not-used
                                       pres    (:,:,:,n),      & ! [OUT] not-used
                                       pott    (:,:,:,n),      & ! [IN]
                                       qtrc    (:,:,:,n,I_QV), & ! [IN]
                                       qc      (:,:,:),        & ! [IN]
                                       temp_sfc(:,:,:,n),      & ! [OUT] not-used
                                       pres_sfc(:,:,:,n),      & ! [IN]
                                       pott_sfc(:,:,:,n),      & ! [IN]
                                       qtrc_sfc(:,:,:,n,I_QV), & ! [IN]
                                       qc_sfc  (:,:,:)         ) ! [IN]

       call COMM_vars8( dens(:,:,:,n), 3 )
       call COMM_wait ( dens(:,:,:,n), 3 )
    end do

    do n = start_step, end_step
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE-1
          momz(k,i,j,n) = velz(k,i,j,n) * ( dens(k+1,i  ,j  ,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do i = 1, IA
          momz(k,i,j,n) = 0.0_RP
       end do
       end do
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          rhot(k,i,j,n) = pott(k,i,j,n) * dens(k,i,j,n)
       end do
       end do
       end do

       do j = 1, JA
       do i = 1, IA-1
       do k = KS, KE
          momx(k,i,j,n) = velx(k,i,j,n) * ( dens(k  ,i+1,j  ,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do k = KS, KE
          momx(k,IA,j,n) = velx(k,i,j,n) * dens(k,i,j,n)
       end do
       end do
       call COMM_vars8( momx(:,:,:,n), 1 )
       call COMM_wait ( momx(:,:,:,n), 1, .false. )

       do j = 1, JA-1
       do i = 1, IA
       do k = KS, KE
          momy(k,i,j,n) = vely(k,i,j,n) * ( dens(k  ,i  ,j+1,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          momy(k,i,JA,n) = vely(k,i,j,n) * dens(k,i,j,n)
       end do
       end do
       call COMM_vars8( momy(:,:,:,n), 1 )
       call COMM_wait ( momy(:,:,:,n), 1, .false. )

    end do

    return
  end subroutine InputAtomNICAM

  !-----------------------------------------------------------------------------
  subroutine InputSurfaceSCALE( &
      tg,                  & ! (out)
      strg,                & ! (out)
      roff,                & ! (out)
      qvef,                & ! (out)
      tw,                  & ! (out)
      lst,                 & ! (out)
      ust,                 & ! (out)
      sst,                 & ! (out)
      albw,                & ! (out)
      albg,                & ! (out)
      z0w,                 & ! (out)
      skint,               & ! (out)
      skinw,               & ! (out)
      snowq,               & ! (out)
      snowt,               & ! (out)
      basename_org,        & ! (in)
      use_file_landwater,  & ! (in)
      init_landwater_ratio ) ! (in)
    use scale_const, only: &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, &
       EPS   => CONST_EPS
    use scale_landuse, only: &
       lsmask_nest => LANDUSE_frac_land,  &
       fact_ocean  => LANDUSE_fact_ocean, &
       fact_land   => LANDUSE_fact_land,  &
       fact_urban  => LANDUSE_fact_urban
    use scale_grid_nest, only: &
       PARENT_KMAX,     &
       PARENT_IMAX,     &
       PARENT_JMAX,     &
       PARENT_LKMAX,    &
       NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y, &
       NEST_TILE_ID
    use mod_land_vars, only: &
       convert_WS2VWC
    implicit none

    real(RP), intent(out) :: tg   (:,:,:)
    real(RP), intent(out) :: strg (:,:,:)
    real(RP), intent(out) :: roff (:,:)
    real(RP), intent(out) :: qvef (:,:)
    real(RP), intent(out) :: tw   (:,:)
    real(RP), intent(out) :: lst  (:,:)
    real(RP), intent(out) :: ust  (:,:)
    real(RP), intent(out) :: sst  (:,:)
    real(RP), intent(out) :: albw (:,:,:)
    real(RP), intent(out) :: albg (:,:,:)
    real(RP), intent(out) :: z0w  (:,:)
    real(RP), intent(out) :: skint(:,:)
    real(RP), intent(out) :: skinw(:,:)
    real(RP), intent(out) :: snowq(:,:)
    real(RP), intent(out) :: snowt(:,:)

    character(LEN=*), intent(in) :: basename_org
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    real(RP),         intent(in) :: init_landwater_ratio ! Ratio of land water to storage is constant

    integer, parameter    :: handle = 1

    ! work
    real(RP), allocatable :: read1D(:)
    real(RP), allocatable :: read2D(:,:)
    real(RP), allocatable :: read3D(:,:,:)

    real(RP), allocatable :: lon_org(:,:,:)
    real(RP), allocatable :: lat_org(:,:,:)
    real(RP), allocatable :: lz_org (:,:,:,:)

    real(RP), allocatable :: tg_org    (:,:,:)
    real(RP), allocatable :: strg_org  (:,:,:)
    real(RP), allocatable :: tw_org    (:,:)
    real(RP), allocatable :: lst_org   (:,:)
    real(RP), allocatable :: ust_org   (:,:)
    real(RP), allocatable :: sst_org   (:,:)
    real(RP), allocatable :: albw_org  (:,:,:)
    real(RP), allocatable :: albg_org  (:,:,:)
    real(RP), allocatable :: z0w_org   (:,:)
    real(RP), allocatable :: skint_org (:,:)
    real(RP), allocatable :: skinw_org (:,:)
    real(RP), allocatable :: snowq_org (:,:)
    real(RP), allocatable :: snowt_org (:,:)
    real(RP), allocatable :: lsmask_org(:,:)

    real(RP), allocatable :: sh2o (:,:,:)

    real(RP), allocatable :: hfact(:,:,:)
    real(RP), allocatable :: vfact(:,:,:,:,:)

    integer, allocatable :: igrd(:,:,:)
    integer, allocatable :: jgrd(:,:,:)
    integer, allocatable :: kgrd(:,:,:,:,:)

    real(RP) :: lcz_3D(LKMAX,IA,JA)

    real(RP)              :: maskval_tg
    real(RP)              :: maskval_strg
    real(RP)              :: maskval_lst
    real(RP), allocatable :: work (:,:,:)

    integer :: KALL, IALL, JALL
    integer :: rank

    integer :: k, i, j, n
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye
    !---------------------------------------------------------------------------

    KALL = PARENT_LKMAX(handle)
    IALL = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    JALL = PARENT_JMAX(handle) * NEST_TILE_NUM_Y

    allocate( read1D(                                           PARENT_LKMAX(handle) ) )
    allocate( read2D( PARENT_IMAX(handle), PARENT_JMAX(handle)                       ) )
    allocate( read3D( PARENT_IMAX(handle), PARENT_JMAX(handle), PARENT_LKMAX(handle) ) )

    allocate( lon_org(       IALL, JALL, 1 ) )
    allocate( lat_org(       IALL, JALL, 1 ) )
    allocate( lz_org ( KALL, IALL, JALL, 1 ) )

    allocate( tg_org    ( KALL, IALL, JALL    ) )
    allocate( strg_org  ( KALL, IALL, JALL    ) )
    allocate( sh2o      ( KALL, IALL, JALL    ) )
    allocate( tw_org    (       IALL, JALL    ) )
    allocate( lst_org   (       IALL, JALL    ) )
    allocate( ust_org   (       IALL, JALL    ) )
    allocate( sst_org   (       IALL, JALL    ) )
    allocate( albw_org  (       IALL, JALL, 2 ) )
    allocate( albg_org  (       IALL, JALL, 2 ) )
    allocate( z0w_org   (       IALL, JALL    ) )
    allocate( skint_org (       IALL, JALL    ) )
    allocate( skinw_org (       IALL, JALL    ) )
    allocate( snowq_org (       IALL, JALL    ) )
    allocate( snowt_org (       IALL, JALL    ) )
    allocate( lsmask_org(       IALL, JALL    ) )
    allocate( work      (       IALL, JALL, 1 ) )

    allocate( hfact(        IA, JA, itp_nh         ) )
    allocate( vfact( LKMAX, IA, JA, itp_nh, itp_nv ) )
    allocate( igrd (        IA, JA, itp_nh         ) )
    allocate( jgrd (        IA, JA, itp_nh         ) )
    allocate( kgrd ( LKMAX, IA, JA, itp_nh, itp_nv ) )

    if( IO_L ) write(IO_FID_LOG,*) ''
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

       call FileRead( read2D(:,:), BASENAME_ORG, "lon", 1, rank )
       lon_org(xs:xe,ys:ye,1) = read2D(:,:) * D2R

       call FileRead( read2D(:,:), BASENAME_ORG, "lat", 1, rank )
       lat_org(xs:xe,ys:ye,1) = read2D(:,:) * D2R

       call FileRead( read1D(:),   BASENAME_ORG, "lz",  1, rank )
       lz_org(:,1,1,1) = read1D(:)

       call FileRead( read3D(:,:,:), BASENAME_ORG, "LAND_TEMP",  1, rank )
       do k = 1, KALL
         tg_org(k,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       if( use_file_landwater )then
        call FileRead( read3D(:,:,:), BASENAME_ORG, "LAND_WATER", 1, rank )
        do k = 1, KALL
          strg_org(k,xs:xe,ys:ye) = read3D(:,:,k)
        end do
       endif

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_TEMP",     1, rank )
       tw_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "LAND_SFC_TEMP",  1, rank )
       lst_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "URBAN_SFC_TEMP", 1, rank )
       ust_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_SFC_TEMP", 1, rank )
       sst_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_ALB_LW",   1, rank )
       albw_org(xs:xe,ys:ye,1) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_ALB_SW",   1, rank )
       albw_org(xs:xe,ys:ye,2) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "LAND_ALB_LW",    1, rank )
       albg_org(xs:xe,ys:ye,1) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "LAND_ALB_SW",    1, rank )
       albg_org(xs:xe,ys:ye,2) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_SFC_Z0M",  1, rank )
       z0w_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "SFC_TEMP",       1, rank )
       skint_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "lsmask",         1, rank )
       lsmask_org(xs:xe,ys:ye) = read2D(:,:)

    end do

    ! fill grid data
    do j = 1, JALL
    do i = 1, IALL
      lz_org(:,i,j,1) = lz_org(:,1,1,1)
    end do
    end do

    do j = 1, JA
    do i = 1, IA
      lcz_3D(:,i,j) = LCZ(:)
    enddo
    enddo

    call latlonz_interpolation_fact( hfact  (:,:,:),     & ! [OUT]
                                     vfact  (:,:,:,:,:), & ! [OUT]
                                     kgrd   (:,:,:,:,:), & ! [OUT]
                                     igrd   (:,:,:),     & ! [OUT]
                                     jgrd   (:,:,:),     & ! [OUT]
                                     lcz_3D (:,:,:),     & ! [IN]
                                     LAT    (:,:),       & ! [IN]
                                     LON    (:,:),       & ! [IN]
                                     lz_org (:,:,:,:),   & ! [IN]
                                     lat_org(:,:,:),     & ! [IN]
                                     lon_org(:,:,:),     & ! [IN]
                                     KALL, IALL, JALL,   & ! [IN]
                                     1,                  & ! [IN]
                                     landgrid=.true.     ) ! [IN]


    ! replace missing value
     maskval_tg   = 298.0_RP    ! mask value 50K => 298K
     maskval_strg = 0.02_RP     ! mask value 0.0 => 0.02
                                ! default value 0.02: set as value of forest at 40% of evapolation rate.
                                ! forest is considered as a typical landuse over Japan area.

    ! Land temp: interpolate over the ocean
     do k = 1, KALL
        work(:,:,1) = tg_org(k,:,:)
        call interp_OceanLand_data(work(:,:,:),lsmask_org, IALL, JALL, 1, landdata=.true.)
        tg_org(k,:,:) = work(:,:,1)
     enddo

    ! interpolation
    do j = 1, JA
    do i = 1, IA
    do k = 1, LKMAX-1
       tg  (k,i,j) = tg_org  (kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                   + tg_org  (kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                   + tg_org  (kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                   + tg_org  (kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                   + tg_org  (kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                   + tg_org  (kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo ! k
       tg  (LKMAX,i,j) = tg  (LKMAX-1,i,j)
    enddo ! i
    enddo ! j
    ! replace values over the ocean
    do k = 1, LKMAX
       call replace_misval( tg(k,:,:), maskval_tg, lsmask_nest)
    enddo

    ! Land water: interpolate over the ocean
    if( use_file_landwater )then
      do k = 1, KALL
         work(:,:,1) = strg_org(k,:,:)
         call interp_OceanLand_data(work(:,:,:), lsmask_org, IALL, JALL, 1, landdata=.true.)
         strg_org(k,:,:) = work(:,:,1)
      enddo
      ! interpolation
      do j = 1, JA
      do i = 1, IA
      do k = 1, LKMAX-1
         strg(k,i,j) = strg_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + strg_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + strg_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + strg_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + strg_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + strg_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,2)
      enddo
         strg(LKMAX,i,j) = strg(LKMAX-1,i,j)
      enddo
      enddo
      ! replace values over the ocean
      do k = 1, LKMAX
       call replace_misval( strg(k,:,:), maskval_strg, lsmask_nest )
      enddo
    else  ! not read from boundary file
      sh2o(:,:,:) = init_landwater_ratio
      ! conversion from water saturation [fraction] to volumetric water content [m3/m3]
      do k = 1, LKMAX
         strg(k,:,:) = convert_WS2VWC( sh2o(k,:,:), critical=.true. )
      end do
    endif

    ! Ocean temp: interpolate over the land
     work(:,:,1) = tw_org(:,:)
     call interp_OceanLand_data(work(:,:,:), lsmask_org, IALL, JALL, 1, landdata=.false.)
     tw_org(:,:) = work(:,:,1)
    ! SST: interpolate over the land
     work(:,:,1) = sst_org(:,:)
     call interp_OceanLand_data(work(:,:,:), lsmask_org, IALL, JALL, 1, landdata=.false.)
     sst_org(:,:) = work(:,:,1)
    ! Surface skin temp: interpolate over the ocean
     work(:,:,1) = lst_org(:,:)
     call interp_OceanLand_data(work(:,:,:), lsmask_org, IALL, JALL, 1, landdata=.true.)
     lst_org(:,:) = work(:,:,1)
    ! Urban surface temp: interpolate over the ocean
     !work(:,:,1) = ust_org(:,:)
     !call interp_OceanLand_data(work(:,:,:), lsmask_org, IALL, JALL, 1, landdata=.true.)
     !ust_org(:,:) = work(:,:,1)
     ust_org(:,:) = lst_org(:,:)

    ! cold start approach
    !skint_org(:,:)      = lst_org(:,:)
    skinw_org(:,:)      = 0.0_RP
    snowq_org(:,:)      = 0.0_RP
    snowt_org(:,:)      = TEM00
    albw_org (:,:,I_LW) = 0.04_RP  ! emissivity of water surface : 0.96
    albw_org (:,:,I_SW) = 0.10_RP
    albg_org (:,:,I_LW) = 0.03_RP  ! emissivity of general ground surface : 0.95-0.98
    albg_org (:,:,I_SW) = 0.22_RP
    z0w_org  (:,:)      = 0.001_RP
    roff(:,:) = 0.0_RP ! not necessary
    qvef(:,:) = 0.0_RP ! not necessary

    do j = 1, JA
    do i = 1, IA
       tw   (i,j) = tw_org   (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + tw_org   (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + tw_org   (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       sst  (i,j) = sst_org  (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + sst_org  (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + sst_org  (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       z0w  (i,j) = z0w_org  (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + z0w_org  (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + z0w_org  (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       lst  (i,j) = lst_org  (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + lst_org  (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + lst_org  (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       ust  (i,j) = ust_org  (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + ust_org  (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + ust_org  (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       !skint(i,j) = skint_org(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
       !           + skint_org(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
       !           + skint_org(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       skinw(i,j) = skinw_org(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + skinw_org(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + skinw_org(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       snowq(i,j) = snowq_org(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + snowq_org(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + snowq_org(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       snowt(i,j) = snowt_org(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + snowt_org(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + snowt_org(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
    enddo
    enddo

    ! replace values over the ocean ####
     do j = 1, JA
     do i = 1, IA
        if( abs(lsmask_nest(i,j)-0.0_RP) < EPS )then ! ocean grid
           lst(i,j)   = sst(i,j)
           ust(i,j)   = sst(i,j)
        endif
           skint(i,j) = fact_ocean(i,j) * sst(i,j) &
                      + fact_land (i,j) * lst(i,j) &
                      + fact_urban(i,j) * ust(i,j)
     enddo
     enddo

    do n = 1, 2 ! 1:LW, 2:SW
    do j = 1, JA
    do i = 1, IA
       albw(i,j,n) = albw_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                   + albw_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                   + albw_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
       albg(i,j,n) = albg_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                   + albg_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                   + albg_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    enddo

    deallocate( read1D )
    deallocate( read2D )
    deallocate( read3D )

    deallocate( lon_org )
    deallocate( lat_org )
    deallocate( lz_org  )

    deallocate( tg_org    )
    deallocate( strg_org  )
    deallocate( sh2o      )
    deallocate( tw_org    )
    deallocate( lst_org   )
    deallocate( ust_org   )
    deallocate( sst_org   )
    deallocate( albw_org  )
    deallocate( albg_org  )
    deallocate( z0w_org   )
    deallocate( skint_org )
    deallocate( skinw_org )
    deallocate( snowq_org )
    deallocate( snowt_org )
    deallocate( lsmask_org )

    deallocate( hfact )
    deallocate( vfact )
    deallocate( igrd  )
    deallocate( jgrd  )
    deallocate( kgrd  )

    return
  end subroutine InputSurfaceSCALE

  !-----------------------------------------------------------------------------
  subroutine InputSurfaceWRF( &
      tg,                   & ! (out)
      strg,                 & ! (out)
      roff,                 & ! (out)
      qvef,                 & ! (out)
      tw,                   & ! (out)
      lst,                  & ! (out)
      ust,                  & ! (out)
      sst,                  & ! (out)
      albw,                 & ! (out)
      albg,                 & ! (out)
      z0w,                  & ! (out)
      skint,                & ! (out)
      skinw,                & ! (out)
      snowq,                & ! (out)
      snowt,                & ! (out)
      basename,             & ! (in)
      dims,                 & ! (in)
      use_file_landwater,   & ! (in)
      init_landwater_ratio, & ! (in)
      mdlid,                & ! (in)
      serial                & ! (in)
      )
    use mod_land_vars, only: &
      convert_WS2VWC
    use scale_landuse, only: &
      frac_land  => LANDUSE_frac_land
    implicit none

    real(RP),         intent(out)  :: tg(:,:,:)
    real(RP),         intent(out)  :: strg(:,:,:)
    real(RP),         intent(out)  :: roff(:,:)
    real(RP),         intent(out)  :: qvef(:,:)
    real(RP),         intent(out)  :: tw(:,:)
    real(RP),         intent(out)  :: lst(:,:)
    real(RP),         intent(out)  :: ust(:,:)
    real(RP),         intent(out)  :: sst(:,:)
    real(RP),         intent(out)  :: albw(:,:,:)
    real(RP),         intent(out)  :: albg(:,:,:)
    real(RP),         intent(out)  :: z0w(:,:)
    real(RP),         intent(out)  :: skint(:,:)
    real(RP),         intent(out)  :: skinw(:,:)
    real(RP),         intent(out)  :: snowq(:,:)
    real(RP),         intent(out)  :: snowt(:,:)
    character(LEN=*), intent( in)  :: basename
    integer,          intent( in)  :: mdlid
    integer,          intent( in)  :: dims(:)
    logical,          intent( in)  :: use_file_landwater   ! use land water data from files
    real(RP),         intent( in)  :: init_landwater_ratio ! Ratio of land water to storage is constant,
    logical,          intent( in)  :: serial

    real(RP), allocatable :: org_3D(:,:,:)
    real(RP), allocatable :: org_4D(:,:,:,:)
    real(RP), allocatable :: lat_ORG(:,:,:)
    real(RP), allocatable :: lon_ORG(:,:,:)
    real(RP), allocatable :: zs_org(:,:,:,:)
    real(RP), allocatable :: dzs_org(:,:)
    real(RP), allocatable :: landmask(:,:)
    real(RP), allocatable :: sh2o(:,:,:)

    real(SP), allocatable :: dummy_2D(:,:)       ! for WRF restart file
    real(SP), allocatable :: dummy_3D(:,:,:)     ! for WRF restart file
    real(SP), allocatable :: dummy_4D(:,:,:,:)   ! for WRF restart file

    real(RP), allocatable :: lcz_3D(:,:,:)
    real(RP), allocatable :: hfact(:,:,:)
    real(RP), allocatable :: vfact(:,:,:,:,:)

    integer, allocatable :: igrd(:,:,:)
    integer, allocatable :: jgrd(:,:,:)
    integer, allocatable :: kgrd(:,:,:,:,:)

    real(RP) :: d2r
    real(RP) :: maskval
    integer  :: n, k, i, j, iq, iqw
    integer  :: tcount = 1

    logical  :: do_read
    logical  :: existence

    intrinsic shape
    !---------------------------------------------------------------------------
    d2r = pi / 180.0_RP

    allocate( lcz_3D(LKMAX,IA,JA) )
    allocate( sh2o(LKMAX,IA,JA) )
    allocate( hfact(IA,JA,itp_nh) )
    allocate( vfact(LKMAX,IA,JA,itp_nh,itp_nv) )
    allocate( igrd(IA,JA,itp_nh) )
    allocate( jgrd(IA,JA,itp_nh) )
    allocate( kgrd(LKMAX,IA,JA,itp_nh,itp_nv) )

    allocate( landmask(dims(2),dims(3)) )
    allocate( org_3D(dims(2),dims(3),tcount) )
    allocate( org_4D(dims(7),dims(2),dims(3),tcount) )
    allocate( lat_org(dims(2),dims(3),tcount) )
    allocate( lon_org(dims(2),dims(3),tcount) )
    allocate( zs_org(dims(7),dims(2),dims(3),tcount) )
    allocate( dzs_org(dims(7),tcount) )
    allocate( dummy_2D(dims(7),tcount) )
    allocate( dummy_3D(dims(2),dims(3),tcount) )
    allocate( dummy_4D(dims(7),dims(2),dims(3),tcount) )

    if( serial ) then
       if( myrank == PRC_master ) then
          do_read = .true.
       else
          do_read = .false.
       endif
    else
       do_read = .true.
    endif

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputWRF-Surface]'

    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "XLAT",    1, tcount, myrank, mdlid, single=.true. )
       lat_org(:,:,:) = dummy_3D(:,:,:) * d2r
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "XLONG",   1, tcount, myrank, mdlid, single=.true. )
       lon_org(:,:,:) = dummy_3D(:,:,:) * d2r
       call ExternalFileRead( dummy_2D(:,:),                               &
                      BASENAME, "ZS",      1, tcount, myrank, mdlid, dims(7), single=.true. )
       do n = 1, tcount
       do j = 1, dims(3)
       do i = 1, dims(2)
          zs_org(:,i,j,n) = dummy_2D(:,n)
       enddo
       enddo
       enddo
       call ExternalFileRead( dummy_2D(:,:),                               &
                      BASENAME, "DZS",      1, tcount, myrank, mdlid, dims(7), single=.true. )
       dzs_org(:,:) = dummy_2D(:,:)
    endif

    if( serial ) then
       call COMM_bcast( lat_org(:,:,:),  dims(2), dims(3), tcount          )
       call COMM_bcast( lon_org(:,:,:),  dims(2), dims(3), tcount          )
       call COMM_bcast( zs_org(:,:,:,:), dims(7), dims(2), dims(3), tcount )
       call COMM_bcast( dzs_org(:,:),    dims(7), tcount                   )
    endif

    ! preparing for interpolation
    do j = 1, JA
    do i = 1, IA
       lcz_3D(:,i,j) = lcz(:)
    enddo
    enddo
    call latlonz_interpolation_fact( hfact,vfact,kgrd,igrd,jgrd,lcz_3D,lat,lon, &
                                      zs_org,lat_org,lon_org,dims(7),dims(2),dims(3),1,landgrid=.true. )

    ! land mask (1:land, 0:water)
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "LANDMASK",  1, tcount, myrank, mdlid, single=.true. )
       org_3D(:,:,:) = dummy_3D(:,:,:)
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    n = 1
    do j = 1, dims(3)
    do i = 1, dims(2)
       landmask(i,j) = org_3D(i,j,n)
    enddo
    enddo

    ! soil temperature [K]
    if( do_read ) then
       call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                      BASENAME, "TSLB",  1, tcount, myrank, mdlid, single=.true., landgrid=.true. )
       org_4D(:,:,:,:) = dummy_4D(:,:,:,:)
    endif
    if( serial ) call COMM_bcast( org_4D(:,:,:,:), dims(7),dims(2),dims(3),tcount )
    ! interpolation over the ocean
    do k = 1, dims(7)
     call interp_OceanLand_data(org_4D(k,:,:,:),landmask,dims(2),dims(3),tcount,landdata=.true.,maskval=maskval)
    enddo
    n = 1
    do j = 1, JA
    do i = 1, IA
    do k = 1, LKMAX-1  ! interpolation
       tg(k,i,j) =  org_4D(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                  + org_4D(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                  + org_4D(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                  + org_4D(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                  + org_4D(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                  + org_4D(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    tg(LKMAX,i,j) = tg(LKMAX-1,i,j)
    enddo
    enddo
    do k = 1, LKMAX
     call replace_misval( tg(k,:,:), maskval, frac_land )
    enddo

    ! soil liquid water [m3 m-3] (no wrfout-default)
    if( use_file_landwater ) then
     call ExternalFileVarExistence( existence, BASENAME, "SH2O", myrank, mdlid, single=.true. )
     if ( existence ) then
       if( do_read ) then
          call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                         BASENAME, "SH2O", 1, tcount, myrank, mdlid, single=.true., landgrid=.true.  )
          org_4D(:,:,:,:) = dummy_4D(:,:,:,:)
       endif
       if( serial ) call COMM_bcast( org_4D(:,:,:,:), dims(7),dims(2),dims(3),tcount )
       ! interpolation over the ocean
       do k = 1, dims(7)
        call interp_OceanLand_data(org_4D(k,:,:,:),landmask,dims(2),dims(3),tcount,landdata=.true.,maskval=maskval)
       enddo
       n = 1
       do j = 1, JA
       do i = 1, IA
       do k = 1, LKMAX-1  ! interpolation
          sh2o(k,i,j) =  org_4D(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                       + org_4D(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                       + org_4D(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                       + org_4D(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                       + org_4D(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                       + org_4D(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
       enddo
       sh2o(LKMAX,i,j) = sh2o(LKMAX-1,i,j)
       enddo
       enddo
       do k = 1, LKMAX
        call replace_misval( sh2o(k,:,:), maskval, frac_land )
       enddo
       ! conversion from water saturation [fraction] to volumetric water content [m3/m3]
       do k = 1, LKMAX
          strg(k,:,:) = convert_WS2VWC( sh2o(k,:,:), critical=.true. )
       end do
     else
       ! default value: set as value of forest at 40% of evapolation rate.
       ! forest is considered as a typical landuse over Japan area.
       strg(:,:,:) = 0.02_DP
     endif
    else  ! not read from boundary file
       sh2o(:,:,:) = init_landwater_ratio
       ! conversion from water saturation [fraction] to volumetric water content [m3/m3]
       do k = 1, LKMAX
          strg(k,:,:) = convert_WS2VWC( sh2o(k,:,:), critical=.true. )
       end do
    endif

    ! surface runoff [mm]
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "SFROFF",  1, tcount, myrank, mdlid, single=.true. )
       do n = 1, tcount
       do j = 1, dims(3)
       do i = 1, dims(2)
          org_3D(i,j,n) = dummy_3D(i,j,n) * 1000.0_DP * dwatr
       enddo
       enddo
       enddo
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    n = 1
    do j = 1, JA
    do i = 1, IA
       roff(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                  + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                  + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ! SEA SURFACE TEMPERATURE [K]
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "SST",  1, tcount, myrank, mdlid, single=.true. )
       org_3D(:,:,:) = dummy_3D(:,:,:)
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    ! interpolate over the land
    call interp_OceanLand_data(org_3D(:,:,:),landmask,dims(2),dims(3),tcount,landdata=.false.)
    n = 1
    do j = 1, JA
    do i = 1, IA
       tw(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    !sst(:,:) = lst(:,:)
    sst(:,:) = tw(:,:)

    ! SURFACE SKIN TEMPERATURE [K]
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "TSK",  1, tcount, myrank, mdlid, single=.true. )
       org_3D(:,:,:) = dummy_3D(:,:,:)
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    ! interpolate over the ocean
    call interp_OceanLand_data(org_3D(:,:,:),landmask,dims(2),dims(3),tcount,landdata=.true.,maskval=maskval)
    n = 1
    do j = 1, JA
    do i = 1, IA
       lst(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                 + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                 + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    call replace_misval( lst(:,:), maskval, frac_land )
    ust(:,:) = lst(:,:)
    qvef(:,:) = 0.0_DP ! currently not used

    ! ALBEDO [-]
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "ALBEDO",  1, tcount, myrank, mdlid, single=.true. )
       org_3D(:,:,:) = dummy_3D(:,:,:)
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    n = 1
    do j = 1, JA
    do i = 1, IA
       albw(i,j,I_SW) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                       + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                       + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    albg(:,:,I_SW) = albw(:,:,I_SW)

    ! SURFACE EMISSIVITY [-]
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "EMISS",  1, tcount, myrank, mdlid, single=.true. )
       do n = 1, tcount
       do j = 1, dims(3)
       do i = 1, dims(2)
          org_3D(i,j,n) = 1.0_DP - dummy_3D(i,j,n)
       enddo
       enddo
       enddo
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    n = 1
    do j = 1, JA
    do i = 1, IA
       albw(i,j,I_LW) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                       + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                       + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    albg(:,:,I_LW) = albw(:,:,I_LW)

    ! TIME-VARYING ROUGHNESS LENGTH [m] (no wrfout-default)
    call ExternalFileVarExistence( existence, BASENAME, "ZNT", myrank, mdlid, single=.true. )
    if ( existence ) then
       if( do_read ) then
          call ExternalFileRead( dummy_3D(:,:,:),                             &
                         BASENAME, "ZNT",  1, tcount, myrank, mdlid, single=.true. )
          org_3D(:,:,:) = dummy_3D(:,:,:)
       endif
       if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
       n = 1
       do j = 1, JA
       do i = 1, IA
          z0w(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                    + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                    + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
       enddo
       enddo
    else
       z0w(:,:) = 0.1_DP
    endif

    !tentative approach for skin
    skint(:,:) = lst(:,:)
    skinw(:,:) = 0.0_DP

    ! SNOW WATER EQUIVALENT [kg m-2] (no wrfout-default)
    call ExternalFileVarExistence( existence, BASENAME, "SNOW", myrank, mdlid, single=.true. )
    if ( existence ) then
       if( do_read ) then
          call ExternalFileRead( dummy_3D(:,:,:),                             &
                         BASENAME, "SNOW",  1, tcount, myrank, mdlid, single=.true. )
          org_3D(:,:,:) = dummy_3D(:,:,:)
       endif
       if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
       n = 1
       do j = 1, JA
       do i = 1, IA
          snowq(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                      + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                      + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
       enddo
       enddo
    else
       snowq(:,:) = 0.0_DP
    endif

    ! AVERAGE SNOW TEMPERATURE [C] (no wrfout-default)
    call ExternalFileVarExistence( existence, BASENAME, "TSNAV", myrank, mdlid, single=.true. )
    if ( existence ) then
       if( do_read ) then
          call ExternalFileRead( dummy_3D(:,:,:),                             &
                         BASENAME, "TSNAV",  1, tcount, myrank, mdlid, single=.true. )
          do n = 1, tcount
          do j = 1, dims(3)
          do i = 1, dims(2)
             org_3D(i,j,n) = dummy_3D(i,j,n) + TEM00
          enddo
          enddo
          enddo
       endif
       if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
       n = 1
       do j = 1, JA
       do i = 1, IA
          snowt(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                      + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                      + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
       enddo
       enddo
    else
       snowt(:,:) = lst(:,:)
    endif

    deallocate( dummy_2D )
    deallocate( dummy_3D )
    deallocate( dummy_4D )
    deallocate( org_3D )
    deallocate( org_4D )
    deallocate( landmask )

    deallocate( lat_org )
    deallocate( lon_org )
    deallocate( zs_org )
    deallocate( dzs_org )

    deallocate( lcz_3D )
    deallocate( hfact )
    deallocate( vfact )
    deallocate( igrd )
    deallocate( jgrd )
    deallocate( kgrd )

    return
  end subroutine InputSurfaceWRF

  !-----------------------------------------------------------------------------
  subroutine InputSurfaceNICAM( &
      tg,                  & ! (out)
      strg,                & ! (out)
      roff,                & ! (out)
      qvef,                & ! (out)
      tw,                  & ! (out)
      lst,                 & ! (out)
      ust,                 & ! (out)
      sst,                 & ! (out)
      albw,                & ! (out)
      albg,                & ! (out)
      z0w,                 & ! (out)
      skint,               & ! (out)
      skinw,               & ! (out)
      snowq,               & ! (out)
      snowt,               & ! (out)
      basename_num,        & ! (in)
      dims,                & ! (in)
      skiplen,             & ! (in)
      use_file_landwater,  & ! (in)
      init_landwater_ratio ) ! (in)
    use scale_const, only: &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, &
       EPS   => CONST_EPS
    use scale_landuse, only: &
       frac_land  => LANDUSE_frac_land
    use mod_land_vars, only: &
      convert_WS2VWC
    implicit none

    real(RP), intent(out) :: tg   (:,:,:)
    real(RP), intent(out) :: strg (:,:,:)
    real(RP), intent(out) :: roff (:,:)
    real(RP), intent(out) :: qvef (:,:)
    real(RP), intent(out) :: tw   (:,:)
    real(RP), intent(out) :: lst  (:,:)
    real(RP), intent(out) :: ust  (:,:)
    real(RP), intent(out) :: sst  (:,:)
    real(RP), intent(out) :: albw (:,:,:)
    real(RP), intent(out) :: albg (:,:,:)
    real(RP), intent(out) :: z0w  (:,:)
    real(RP), intent(out) :: skint(:,:)
    real(RP), intent(out) :: skinw(:,:)
    real(RP), intent(out) :: snowq(:,:)
    real(RP), intent(out) :: snowt(:,:)

    character(LEN=*), intent(in) :: basename_num
    integer,          intent(in) :: dims(:)
    integer,          intent(in) :: skiplen
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    real(RP),         intent(in) :: init_landwater_ratio ! Ratio of land water to storage is constant,
                                                         !              if use_file_landwater is ".false."

    ! ----------------------------------------------------------------

    ! work
    real(SP), allocatable :: read1DX(:)
    real(SP), allocatable :: read1DY(:)
    real(SP), allocatable :: read1DZ(:)
    real(SP), allocatable :: read3DT(:,:,:,:)
    real(SP), allocatable :: read3DS(:,:,:,:)
    real(SP), allocatable :: read4D (:,:,:,:)

    real(RP), allocatable :: lon_org(:,:,:)
    real(RP), allocatable :: lat_org(:,:,:)
    real(RP), allocatable :: lz_org (:,:,:,:)

    real(RP), allocatable :: landmask (:,:)
    real(RP), allocatable :: tg_org   (:,:,:)
    real(RP), allocatable :: strg_org (:,:,:)
    real(RP), allocatable :: tw_org   (:,:)
    real(RP), allocatable :: lst_org  (:,:)
    real(RP), allocatable :: ust_org  (:,:)
    real(RP), allocatable :: sst_org  (:,:)
    real(RP), allocatable :: ice_org  (:,:)
    real(RP), allocatable :: albw_org (:,:,:)
    real(RP), allocatable :: albg_org (:,:,:)
    real(RP), allocatable :: z0w_org  (:,:)
    real(RP), allocatable :: skint_org(:,:)
    real(RP), allocatable :: skinw_org(:,:)
    real(RP), allocatable :: snowq_org(:,:)
    real(RP), allocatable :: snowt_org(:,:)
    real(RP), allocatable :: sh2o(:,:,:)

    real(RP), allocatable :: hfact(:,:,:)
    real(RP), allocatable :: vfact(:,:,:,:,:)

    integer,  allocatable :: igrd(:,:,:)
    integer,  allocatable :: jgrd(:,:,:)
    integer,  allocatable :: kgrd(:,:,:,:,:)

    real(RP) :: lcz_3D(LKMAX,IA,JA)

    real(RP)              :: maskval_tg
    real(RP)              :: maskval_strg
    real(RP)              :: maskval_lst
    real(RP), parameter   :: sst_missval=273.1506_RP
    real(RP), allocatable :: work (:,:,:)

    integer :: start_step
    integer :: end_step
    integer :: k, i, j, n, nt

    logical :: do_read

    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    ! read data for initial condition
    start_step = 1
    end_step   = 1
    nt = end_step - start_step + 1

    if( myrank == PRC_master ) then
       do_read = .true.
    else
       do_read = .false.
    endif

    if( do_read ) then
       allocate( read1DX( dims(1)                        ) )
       allocate( read1DY( dims(2)                        ) )
       allocate( read1DZ( dims(7)                        ) )
       allocate( read3DT( 1,       dims(1), dims(2), nt  ) )
       allocate( read3DS( 1,       dims(1), dims(2), nt  ) )
       allocate( read4D ( dims(7), dims(1), dims(2), nt  ) )
    endif

    allocate( lon_org(          dims(1), dims(2), nt  ) )
    allocate( lat_org(          dims(1), dims(2), nt  ) )
    allocate( lz_org ( dims(7), dims(1), dims(2), nt  ) )

    allocate( work     (          dims(1), dims(2), 1 ) )
    allocate( landmask (          dims(1), dims(2)    ) )
    allocate( tg_org   ( dims(7), dims(1), dims(2)    ) )
    allocate( strg_org ( dims(7), dims(1), dims(2)    ) )
    allocate( tw_org   (          dims(1), dims(2)    ) )
    allocate( lst_org  (          dims(1), dims(2)    ) )
    allocate( ust_org  (          dims(1), dims(2)    ) )
    allocate( sst_org  (          dims(1), dims(2)    ) )
    allocate( ice_org  (          dims(1), dims(2)    ) )
    allocate( albw_org (          dims(1), dims(2), 2 ) )
    allocate( albg_org (          dims(1), dims(2), 2 ) )
    allocate( z0w_org  (          dims(1), dims(2)    ) )
    allocate( skint_org(          dims(1), dims(2)    ) )
    allocate( skinw_org(          dims(1), dims(2)    ) )
    allocate( snowq_org(          dims(1), dims(2)    ) )
    allocate( snowt_org(          dims(1), dims(2)    ) )

    allocate( hfact(        IA, JA, itp_nh         ) )
    allocate( vfact( LKMAX, IA, JA, itp_nh, itp_nv ) )
    allocate( igrd (        IA, JA, itp_nh         ) )
    allocate( jgrd (        IA, JA, itp_nh         ) )
    allocate( kgrd ( LKMAX, IA, JA, itp_nh, itp_nv ) )

    allocate( sh2o (LKMAX,IA,JA) )

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputNICAM-Surface]'

    if( do_read ) then
       basename = "la_tg"//trim(basename_num)
       call FileRead( read1DX(:), trim(basename), "lon", start_step, 1, single=.true. )
       do j = 1, dims(2)
          lon_org (:,j,start_step)  = read1DX(:) * D2R
       enddo

       call FileRead( read1DY(:), trim(basename), "lat", start_step, 1, single=.true. )
       do i = 1, dims(1)
          lat_org (i,:,start_step)  = read1DY(:) * D2R
       enddo

       call FileRead( read1DZ(:), trim(basename), "lev", start_step, 1, single=.true. )
       do j = 1, dims(2)
       do i = 1, dims(1)
          lz_org(:,i,j,start_step) = read1DZ(:)
       enddo
       enddo
    endif
    call COMM_bcast( lat_org (:,:,:),            dims(1), dims(2), nt )
    call COMM_bcast( lon_org (:,:,:),            dims(1), dims(2), nt )
    call COMM_bcast( lz_org  (:,:,:,:), dims(7), dims(1), dims(2), nt )

    do j = 1, JA
    do i = 1, IA
      lcz_3D(:,i,j) = LCZ(:)
    enddo
    enddo

    call latlonz_interpolation_fact( hfact  (:,:,:),            & ! [OUT]
                                     vfact  (:,:,:,:,:),        & ! [OUT]
                                     kgrd   (:,:,:,:,:),        & ! [OUT]
                                     igrd   (:,:,:),            & ! [OUT]
                                     jgrd   (:,:,:),            & ! [OUT]
                                     lcz_3D (:,:,:),            & ! [IN]
                                     LAT    (:,:),              & ! [IN]
                                     LON    (:,:),              & ! [IN]
                                     lz_org (:,:,:,:),          & ! [IN]
                                     lat_org(:,:,:),            & ! [IN]
                                     lon_org(:,:,:),            & ! [IN]
                                     dims(7), dims(1), dims(2), & ! [IN]
                                     1,                         & ! [IN]
                                     landgrid=.true.            ) ! [IN]
    deallocate( lon_org )
    deallocate( lat_org )
    deallocate( lz_org  )


    if( do_read ) then
       !
       basename = "lsmask"//trim(basename_num)
       call ExternalFileRead( read3DS(:,:,:,:), &
                              trim(basename),   &
                              "lsmask",         &
                              start_step,       &
                              end_step,         &
                              myrank,           &
                              iNICAM,           &
                              single=.true.     )
       landmask(:,:) = real( read3DS(1,:,:,1), kind=RP )

       ! [scale-offset]
       basename = "la_tg"//trim(basename_num)
       call ExternalFileReadOffset( read4D(:,:,:,:),  &
                                    trim(basename),   &
                                    "la_tg",          &
                                    start_step,       &
                                    end_step,         &
                                    myrank,           &
                                    iNICAM,           &
                                    single=.true.     )
       tg_org(:,:,:) = real( read4D(:,:,:,1), kind=RP )

       ! [scale-offset]
       if( use_file_landwater ) then
          basename = "la_wg"//trim(basename_num)
          call ExternalFileReadOffset( read4D(:,:,:,:),  &
                                       trim(basename),   &
                                       "la_wg",          &
                                       start_step,       &
                                       end_step,         &
                                       myrank,           &
                                       iNICAM,           &
                                       single=.true.     )
          strg_org(:,:,:) = real( read4D(:,:,:,1), kind=RP )
       endif

       ! 6hourly data, first "skiplen" steps are skipped
       start_step = skiplen + 1
       end_step   = skiplen + 1
       basename = "ss_tem_sfc"//trim(basename_num)
       call ExternalFileRead( read3DS(:,:,:,:), &
                              trim(basename),   &
                              "ss_tem_sfc",     &
                              start_step,       &
                              end_step,         &
                              myrank,           &
                              iNICAM,           &
                              single=.true.     )
       lst_org(:,:) = real( read3DS(1,:,:,1), kind=RP )

       ! [scale-offset]
       start_step = 1
       end_step   = 1
       basename = "oa_sst"//trim(basename_num)
       call ExternalFileReadOffset( read3DT(:,:,:,:), &
                                    trim(basename),   &
                                    "oa_sst",         &
                                    start_step,       &
                                    end_step,         &
                                    myrank,           &
                                    iNICAM,           &
                                    single=.true.     )
       sst_org(:,:) = real( read3DT(1,:,:,1), kind=RP )

       !basename = "oa_ice"//trim(basename_num)
       !call ExternalFileRead( read3DS(:,:,:,:), trim(basename), "oa_ice", start_step, end_step, myrank, iNICAM, single=.true. )
       !ice_org(:,:) = real( read3DS(1,:,:,1), kind=RP )

       deallocate( read1DX )
       deallocate( read1DY )
       deallocate( read1DZ )
       deallocate( read3DT )
       deallocate( read3DS )
       deallocate( read4D  )
    endif

    call COMM_bcast( landmask(:,:),              dims(1), dims(2)     )
    call COMM_bcast( lst_org (:,:),              dims(1), dims(2)     )
    call COMM_bcast( sst_org (:,:),              dims(1), dims(2)     )
    !call COMM_bcast( ice_org (:,:),              dims(1), dims(2)     )
    call COMM_bcast( tg_org  (:,:,:),   dims(7), dims(1), dims(2)     )
    if( use_file_landwater ) then
     call COMM_bcast( strg_org(:,:,:),   dims(7), dims(1), dims(2)     )
    endif

    ! replace missing value
     maskval_tg   = 298.0_RP    ! mask value 50K => 298K
     maskval_strg = 0.02_RP     ! mask value 0.0 => 0.02
                                ! default value 0.02: set as value of forest at 40% of evapolation rate.
                                ! forest is considered as a typical landuse over Japan area.

    ! Land temp: interpolate over the ocean
     do k=1,dims(7)
     do j=1,dims(2)
     do i=1,dims(1)
        if ( abs(tg_org(k,i,j)-50.0D0) < EPS ) then
           tg_org(k,i,j) = maskval_tg
        endif
     enddo
     enddo
     enddo

     do k=1,dims(7)
        work(:,:,1) = tg_org(k,:,:)
        call interp_OceanLand_data(work(:,:,:),landmask,dims(1),dims(2), 1,landdata=.true.)
        tg_org(k,:,:) = work(:,:,1)
     enddo
     ! interpolation
     do j = 1, JA
     do i = 1, IA
     do k = 1, LKMAX-1
        tg  (k,i,j) = tg_org  (kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                    + tg_org  (kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                    + tg_org  (kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                    + tg_org  (kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                    + tg_org  (kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                    + tg_org  (kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,2)
     enddo
        tg  (LKMAX,i,j) = tg  (LKMAX-1,i,j)
     enddo
     enddo
     deallocate( tg_org )
     ! replace values over the ocean
     do k = 1, LKMAX
      call replace_misval( tg(k,:,:),   maskval_tg,   frac_land )
     enddo

    ! Land water: interpolate over the ocean
    if( use_file_landwater ) then
      do j=1,dims(2)
      do i=1,dims(1)
      do k=1,dims(7)
         if ( abs(strg_org(k,i,j)-0.0D0) < EPS ) then
            strg_org(k,i,j) = maskval_strg
         endif
      enddo
      enddo
      enddo
      do k=1,dims(7)
         work(:,:,1) = strg_org(k,:,:)
         call interp_OceanLand_data(work(:,:,:),landmask,dims(1),dims(2), 1,landdata=.true.)
         strg_org(k,:,:) = work(:,:,1)
      enddo
      ! interpolation
      do j = 1, JA
      do i = 1, IA
      do k = 1, LKMAX-1
         strg(k,i,j) = strg_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                     + strg_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                     + strg_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                     + strg_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                     + strg_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                     + strg_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,2)
      enddo
         strg(LKMAX,i,j) = strg(LKMAX-1,i,j)
      enddo
      enddo
      ! replace values over the ocean
      do k = 1, LKMAX
       call replace_misval( strg(k,:,:), maskval_strg, frac_land )
      enddo
    else
      sh2o(:,:,:) = init_landwater_ratio
      ! conversion from water saturation [fraction] to volumetric water content [m3/m3]
      do k = 1, LKMAX
         strg(k,:,:) = convert_WS2VWC( sh2o(k,:,:), critical=.true. )
      end do
    endif
    deallocate( strg_org )

    ! Surface skin temp: interpolate over the ocean
     work(:,:,1) = lst_org(:,:)
     call interp_OceanLand_data(work(:,:,:),landmask,dims(1),dims(2), 1,landdata=.true.)
     lst_org(:,:) = work(:,:,1)

    ! SST: retrieve SST data around coast
     do j=1,dims(2)
     do i=1,dims(1)
        if( abs(landmask(i,j)-1.0_RP) < EPS )then  ! not land data
           cycle
        else
           sst_org(i,j)=(sst_org(i,j)-sst_missval*landmask(i,j))/(1.0_RP-landmask(i,j))
        endif
     enddo
     enddo
    ! interpolate over the land
     work(:,:,1) = sst_org(:,:)
     call interp_OceanLand_data(work(:,:,:),landmask,dims(1),dims(2), 1,landdata=.false.)
     sst_org(:,:) = work(:,:,1)

    ! cold start approach
    skint_org(:,:)      = lst_org(:,:)
    skinw_org(:,:)      = 0.0_RP
    snowq_org(:,:)      = 0.0_RP
    snowt_org(:,:)      = TEM00
    tw_org   (:,:)      = sst_org(:,:)
    ust_org  (:,:)      = lst_org(:,:)
    albw_org (:,:,I_LW) = 0.04_RP  ! emissivity of water surface : 0.96
    albw_org (:,:,I_SW) = 0.10_RP
    albg_org (:,:,I_LW) = 0.03_RP  ! emissivity of general ground surface : 0.95-0.98
    albg_org (:,:,I_SW) = 0.22_RP
    z0w_org  (:,:)      = 0.001_RP

    roff(:,:) = 0.0_RP ! not necessary
    qvef(:,:) = 0.0_RP ! not necessary

    do j = 1, JA
    do i = 1, IA
       tw   (i,j) = tw_org   (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + tw_org   (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + tw_org   (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       sst  (i,j) = sst_org  (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + sst_org  (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + sst_org  (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       z0w  (i,j) = z0w_org  (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + z0w_org  (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + z0w_org  (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       lst  (i,j) = lst_org  (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + lst_org  (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + lst_org  (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       ust  (i,j) = ust_org  (igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + ust_org  (igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + ust_org  (igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       skint(i,j) = skint_org(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + skint_org(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + skint_org(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       skinw(i,j) = skinw_org(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + skinw_org(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + skinw_org(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       snowq(i,j) = snowq_org(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + snowq_org(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + snowq_org(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
       snowt(i,j) = snowt_org(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) &
                  + snowt_org(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) &
                  + snowt_org(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
    enddo
    enddo

    ! replace values over the ocean
     do j = 1, JA
     do i = 1, IA
        if( abs(frac_land(i,j)-0.0_RP) < EPS )then ! ocean grid
           lst(i,j)   = sst(i,j)
           ust(i,j)   = sst(i,j)
           skint(i,j) = sst(i,j)
        endif
     enddo
     enddo


    do n = 1, 2 ! 1:LW, 2:SW
    do j = 1, JA
    do i = 1, IA
       albw(i,j,n) = albw_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                   + albw_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                   + albw_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
       albg(i,j,n) = albg_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                   + albg_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                   + albg_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    enddo

    deallocate( tw_org    )
    deallocate( lst_org   )
    deallocate( ust_org   )
    deallocate( sst_org   )
    deallocate( ice_org   )
    deallocate( albw_org  )
    deallocate( albg_org  )
    deallocate( z0w_org   )
    deallocate( skint_org )
    deallocate( skinw_org )
    deallocate( snowq_org )
    deallocate( snowt_org )

    deallocate( sh2o )
    deallocate( hfact )
    deallocate( vfact )
    deallocate( igrd  )
    deallocate( jgrd  )
    deallocate( kgrd  )

    return
  end subroutine InputSurfaceNICAM

  !-----------------------------------------------------------------------------
  subroutine latlonz_interpolation_fact( &
      hfact,      & ! (out)
      vfact,      & ! (out)
      kgrd,       & ! (out)
      igrd,       & ! (out)
      jgrd,       & ! (out)
      myhgt,      & ! (in)
      mylat,      & ! (in)
      mylon,      & ! (in)
      inhgt,      & ! (in)
      inlat,      & ! (in)
      inlon,      & ! (in)
      nz,         & ! (in)
      nx,         & ! (in)
      ny,         & ! (in)
      step,       & ! (in)
      landgrid    ) ! (in)
    implicit none

    real(RP), intent(out) :: hfact(:,:,:)
    real(RP), intent(out) :: vfact(:,:,:,:,:)
    integer,  intent(out) :: kgrd (:,:,:,:,:)
    integer,  intent(out) :: igrd (:,:,:)
    integer,  intent(out) :: jgrd (:,:,:)

    real(RP), intent(in)  :: myhgt(:,:,:)
    real(RP), intent(in)  :: mylat(:,:)
    real(RP), intent(in)  :: mylon(:,:)
    real(RP), intent(in)  :: inhgt(:,:,:,:)
    real(RP), intent(in)  :: inlat(:,:,:)
    real(RP), intent(in)  :: inlon(:,:,:)
    integer,  intent(in)  :: nz
    integer,  intent(in)  :: nx
    integer,  intent(in)  :: ny
    integer,  intent(in)  :: step

    logical,  intent(in), optional :: landgrid

    real(RP) :: distance
    real(RP) :: denom
    real(RP) :: dist(itp_nh)
    integer :: i, j, k, ii, jj, kk
    integer :: idx
    integer :: istart, iend, iinc, blk_i
    integer :: jstart, jend, jinc, blk_j
    integer :: kstart, kend
    logical :: lndgrd
    !---------------------------------------------------------------------------

    lndgrd = .false.
    if ( present(landgrid) .and. landgrid ) then
       lndgrd = .true.
    endif

    hfact(:,:,:) = 0.0_RP
    vfact(:,:,:,:,:) = 0.0_RP

    do j = 1, JA
    do i = 1, IA
       ! nearest block search
       iinc = (nx + 1) / interp_search_divnum
       jinc = (ny + 1) / interp_search_divnum
       dist(1) = large_number_one
       jj = 1 + (jinc/2)
       do while (jj <= ny)
          ii = 1 + (iinc/2)
          do while (ii <= nx)
             distance = haversine( mylat(i,j),mylon(i,j),inlat(ii,jj,step),inlon(ii,jj,step) )

             if( distance < dist(1) )then
                dist(1) = distance
                blk_i = ii
                blk_j = jj
             endif
             ii = ii + iinc
          enddo
          jj = jj + jinc
       enddo
       istart = blk_i - (iinc/2) - 1
       if( istart < 1 ) istart = 1
       iend   = blk_i + (iinc/2) + 1
       if( iend  > nx ) iend   = nx
       jstart = blk_j - (jinc/2) - 1
       if( jstart < 1 ) jstart = 1
       jend   = blk_j + (jinc/2) + 1
       if( jend  > ny ) jend   = ny

       ! main search
       dist(1) = large_number_three
       dist(2) = large_number_two
       dist(3) = large_number_one
       do jj = jstart, jend
       do ii = istart, iend
          distance = haversine( mylat(i,j),mylon(i,j),inlat(ii,jj,step),inlon(ii,jj,step) )
          if ( distance <= dist(1) ) then
             dist(3) = dist(2);     igrd(i,j,3) = igrd(i,j,2);  jgrd(i,j,3) = jgrd(i,j,2)
             dist(2) = dist(1);     igrd(i,j,2) = igrd(i,j,1);  jgrd(i,j,2) = jgrd(i,j,1)
             dist(1) = distance;    igrd(i,j,1) = ii;           jgrd(i,j,1) = jj
          elseif ( dist(1) < distance .and. distance <= dist(2) ) then
             dist(3) = dist(2);     igrd(i,j,3) = igrd(i,j,2);  jgrd(i,j,3) = jgrd(i,j,2)
             dist(2) = distance;    igrd(i,j,2) = ii;           jgrd(i,j,2) = jj
          elseif ( dist(2) < distance .and. distance <= dist(3) ) then
             dist(3) = distance;    igrd(i,j,3) = ii;           jgrd(i,j,3) = jj
          endif
       enddo
       enddo
       if( dist(1)==0.0_RP )then
          hfact(i,j,1) = 1.0_RP
          hfact(i,j,2) = 0.0_RP
          hfact(i,j,3) = 0.0_RP
       else
          denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) + (1.0_RP/dist(3)) )
          hfact(i,j,1) = ( 1.0_RP/dist(1) ) * denom
          hfact(i,j,2) = ( 1.0_RP/dist(2) ) * denom
          hfact(i,j,3) = ( 1.0_RP/dist(3) ) * denom
       endif

       if ( lndgrd ) then
          kstart = 1;     kend = LKMAX
       else
          kstart = KS-1;  kend = KE+1
       endif
       do idx = 1, itp_nh
          ii = igrd(i,j,idx)
          jj = jgrd(i,j,idx)
          do k = kstart, kend
             dist(1) = large_number_two
             dist(2) = large_number_one
             do kk = 1, nz
                distance = abs( myhgt(k,i,j) - inhgt(kk,ii,jj,step) )
                if ( distance <= dist(1) ) then
                   dist(2) = dist(1);     kgrd(k,i,j,idx,2) = kgrd(k,i,j,idx,1)
                   dist(1) = distance;    kgrd(k,i,j,idx,1) = kk
                elseif ( dist(1) < distance .and. distance <= dist(2) ) then
                   dist(2) = distance;    kgrd(k,i,j,idx,2) = kk
                endif
             enddo
             if( dist(1)==0.0_RP )then
                vfact(k,i,j,idx,1) = 1.0_RP
                vfact(k,i,j,idx,2) = 0.0_RP
             else
                denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) )
                vfact(k,i,j,idx,1) = ( 1.0_RP/dist(1) ) * denom
                vfact(k,i,j,idx,2) = ( 1.0_RP/dist(2) ) * denom
             endif
          enddo
       enddo
    enddo
    enddo

    return
  end subroutine latlonz_interpolation_fact

  !-----------------------------------------------------------------------------
  subroutine interp_OceanLand_data( &
      data,      & ! (inout)
      lsmask,    & ! (in)
      nx,        & ! (in)
      ny,        & ! (in)
      tcount,    & ! (in)
      landdata,  & ! (in)
      maskval    & ! (in)
      )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none
    real(RP), intent(inout)         :: data(:,:,:)
    real(RP), intent(in)            :: lsmask(:,:)
    integer,  intent(in)            :: nx
    integer,  intent(in)            :: ny
    integer,  intent(in)            :: tcount
    logical,  intent(in)            :: landdata   ! .true. => land data , .false. => ocean data
    real(RP), intent(out), optional :: maskval

    integer                 :: untarget_mask
    integer, allocatable    :: imask(:,:),imaskr(:,:)
    real(RP),allocatable    :: newdata(:,:)
    real(RP)                :: flag(8)
    real(RP)                :: ydata(8)
    real(RP)                :: dd

    integer :: i, j, k, ii, jj, kk ,n

    !---------------------------------------------------------------------------
    allocate( imask  (nx,ny) )
    allocate( imaskr (nx,ny) )
    allocate( newdata(nx,ny) )
    newdata = 0.0_RP
    if( present(maskval) ) maskval = 999.99_RP

    ! search target cell for interpolation
     do j = 1, ny
     do i = 1, nx
        if( abs(lsmask(i,j)-1.0_RP) < EPS )then
           imask(i,j) = 1  ! land grid
        else
           imask(i,j) = 0  ! ocean grid
        endif
     enddo
     enddo
     if ( landdata ) then  ! interpolation for land data
       untarget_mask = 1
     else                  ! interpolation for ocean data
       untarget_mask = 0
     endif

    ! start interpolation
    do n = 1, tcount

    imaskr = imask
    do kk = 1, 20
    do j  = 1, ny
    do i  = 1, nx
       if ( imask(i,j) == untarget_mask ) then  ! not missing value
           newdata(i,j) = data(i,j,n)
           cycle
       else

         if ( present(maskval).and.(abs(maskval-999.99_RP)<EPS) ) then
            if (abs(lsmask(i,j)-0.0_RP) < EPS) maskval = data(i,j,n)
         endif

        !--------------------------------------
        ! check data of neighbor grid
        !
        !           flag(i,j,8)
        !
        !     +---+---+---+---+---+---+
        !     |       |       |       |
        !     +   6   +   7   +   8   +
        !     |       |       |       |
        !     +---+---+---+---+---+---+
        !     |       |       |       |
        !     +   4   + (i,j) +   5   +
        !     |       |       |       |
        !     +---+---+---+---+---+---+
        !     |       |       |       |
        !     +   1   +   2   +   3   +
        !     |       |       |       |
        !     +---+---+---+---+---+---+
        !---------------------------------------
        flag    = 0.0_RP
        ydata   = 0.0_RP
        if((i==1).and.(j==1))then
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
           if(imask(i+1,j+1)==untarget_mask) flag(8)=1.0_RP
        else if((i==nx).and.(j==1))then
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i-1,j+1)==untarget_mask) flag(6)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
        else if((i==1).and.(j==ny))then
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i+1,j-1)==untarget_mask) flag(3)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
        else if((i==nx).and.(j==ny))then
           if(imask(i-1,j-1)==untarget_mask) flag(1)=1.0_RP
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
        else if(i==1)then
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i+1,j-1)==untarget_mask) flag(3)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
           if(imask(i+1,j+1)==untarget_mask) flag(8)=1.0_RP
        else if(i==nx)then
           if(imask(i-1,j-1)==untarget_mask) flag(1)=1.0_RP
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i-1,j+1)==untarget_mask) flag(6)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
        else if(j==1)then
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
           if(imask(i-1,j+1)==untarget_mask) flag(6)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
           if(imask(i+1,j+1)==untarget_mask) flag(8)=1.0_RP
        else if(j==ny)then
           if(imask(i-1,j-1)==untarget_mask) flag(1)=1.0_RP
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i+1,j-1)==untarget_mask) flag(3)=1.0_RP
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
        else
           if(imask(i-1,j-1)==untarget_mask) flag(1)=1.0_RP
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i+1,j-1)==untarget_mask) flag(3)=1.0_RP
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
           if(imask(i-1,j+1)==untarget_mask) flag(6)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
           if(imask(i+1,j+1)==untarget_mask) flag(8)=1.0_RP
        endif

        if(int(sum(flag(:)))>=3)then  ! coast grid : interpolate
           dd = 0.0_RP
           newdata(i,j) = 0.0_RP
           if(int(flag(1))==1) newdata(i,j) = newdata(i,j) + data(i-1,j-1,n)
           if(int(flag(2))==1) newdata(i,j) = newdata(i,j) + data(i,j-1,n)
           if(int(flag(3))==1) newdata(i,j) = newdata(i,j) + data(i+1,j-1,n)
           if(int(flag(4))==1) newdata(i,j) = newdata(i,j) + data(i-1,j,n)
           if(int(flag(5))==1) newdata(i,j) = newdata(i,j) + data(i+1,j,n)
           if(int(flag(6))==1) newdata(i,j) = newdata(i,j) + data(i-1,j+1,n)
           if(int(flag(7))==1) newdata(i,j) = newdata(i,j) + data(i,j+1,n)
           if(int(flag(8))==1) newdata(i,j) = newdata(i,j) + data(i+1,j+1,n)

           dd = sum(flag(:))
           newdata(i,j) = newdata(i,j) / dd

           imaskr(i,j) = untarget_mask
        else
          newdata(i,j) = data(i,j,n)
        endif

       endif ! sea/land

    enddo
    enddo

     imask(:,:) = imaskr(:,:)
     data(:,:,n)  = newdata(:,:)
    enddo ! kk

    enddo ! n

    deallocate( imask   )
    deallocate( imaskr  )
    deallocate( newdata )

    return
  end subroutine interp_OceanLand_data

  !-----------------------------------------------------------------------------
  subroutine replace_misval( orgdata, maskval, frac_land )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none
    real(RP), intent(inout) :: orgdata(:,:)
    real(RP), intent(in)    :: maskval
    real(RP), intent(in)    :: frac_land(:,:)
    integer                 :: i, j

    do j = 1, JA
    do i = 1, IA
       if( abs(frac_land(i,j)-0.0_RP) < EPS )then ! ocean grid
          orgdata(i,j) = maskval
       endif
    enddo
    enddo

  end subroutine replace_misval

  !-----------------------------------------------------------------------------
  subroutine diagnose_number_concentration( &
      qvars      & ! (inout)
      )
    implicit none
    real(RP), intent(inout) :: qvars(:,:,:,:,:)

    real(RP), parameter :: Dc   =  20.D-6  ! typical particle diameter for cloud  [m]
    real(RP), parameter :: Dr   = 200.D-6  ! typical particle diameter for rain   [m]
    real(RP), parameter :: Di   =  80.D-6  ! typical particle diameter for ice    [m]
    real(RP), parameter :: Ds   =  80.D-6  ! typical particle diameter for snow   [m]
    real(RP), parameter :: Dg   = 200.D-6  ! typical particle diameter for grapel [m]
    real(RP), parameter :: RHOw = 1000.D0  ! typical density for warm particles   [kg/m3]
    real(RP), parameter :: RHOf =  100.D0  ! typical density for frozen particles [kg/m3]
    real(RP), parameter :: RHOg =  400.D0  ! typical density for grapel particles [kg/m3]
    real(RP), parameter :: b    =  3.D0    ! assume spherical form

    real(RP) :: piov6
    integer :: i, j, k, t
    !---------------------------------------------------------------------------

    piov6 = pi / 6.0_RP

#ifndef DRY
    qvars(:,:,:,:,I_NC) = qvars(:,:,:,:,I_QC) / ( (piov6*RHOw) * Dc**b )
    qvars(:,:,:,:,I_NR) = qvars(:,:,:,:,I_QR) / ( (piov6*RHOw) * Dr**b )
    qvars(:,:,:,:,I_NI) = qvars(:,:,:,:,I_QI) / ( (piov6*RHOf) * Di**b )
    qvars(:,:,:,:,I_NS) = qvars(:,:,:,:,I_QS) / ( (piov6*RHOf) * Ds**b )
    qvars(:,:,:,:,I_NG) = qvars(:,:,:,:,I_QG) / ( (piov6*RHOg) * Dg**b )
#endif

    return
  end subroutine diagnose_number_concentration

  !-----------------------------------------------------------------------------
  subroutine check_domain_compatibility( &
      lon_org,     & ! (in)
      lat_org,     & ! (in)
      lev_org,     & ! (in)
      lon_loc,     & ! (in)
      lat_loc,     & ! (in)
      lev_loc,     & ! (in)
      skip_x,      & ! (in)
      skip_y,      & ! (in)
      skip_z       ) ! (in)
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none
    real(RP), intent(in) :: lon_org(:,:)
    real(RP), intent(in) :: lat_org(:,:)
    real(RP), intent(in) :: lev_org(:,:,:)
    real(RP), intent(in) :: lon_loc(:,:)
    real(RP), intent(in) :: lat_loc(:,:)
    real(RP), intent(in) :: lev_loc(:,:,:)
    logical,  intent(in), optional :: skip_x
    logical,  intent(in), optional :: skip_y
    logical,  intent(in), optional :: skip_z

    real(RP) :: max_ref, min_ref
    real(RP) :: max_loc, min_loc

    logical :: do_xdirec
    logical :: do_ydirec
    logical :: do_zdirec
    !---------------------------------------------------------------------------

    do_xdirec = .true.
    if ( present(skip_x) .and. skip_x ) then
       do_xdirec = .false.
    endif

    do_ydirec = .true.
    if ( present(skip_y) .and. skip_y ) then
       do_ydirec = .false.
    endif

    do_zdirec = .true.
    if ( present(skip_z) .and. skip_z ) then
       do_zdirec = .false.
    endif

    if ( do_xdirec ) then
       max_ref = maxval( lon_org(:,:) / D2R )
       min_ref = minval( lon_org(:,:) / D2R )
       max_loc = maxval( lon_loc(:,:) / D2R )
       min_loc = minval( lon_loc(:,:) / D2R )

       if ( max_ref < max_loc .or. min_ref > min_loc ) then
          write(*,*) 'xxx ERROR: REQUESTED DOMAIN IS TOO MUCH BROAD'
          write(*,*) 'xxx -- LONGITUDINAL direction over the limit'
          write(*,*) 'xxx -- reference max: ', max_ref
          write(*,*) 'xxx -- reference min: ', min_ref
          write(*,*) 'xxx --     local max: ', max_loc
          write(*,*) 'xxx --     local min: ', min_loc
          call PRC_MPIstop
       endif
    endif

    if ( do_ydirec ) then
       max_ref = maxval( lat_org(:,:) / D2R )
       min_ref = minval( lat_org(:,:) / D2R )
       max_loc = maxval( lat_loc(:,:) / D2R )
       min_loc = minval( lat_loc(:,:) / D2R )

       if ( max_ref < max_loc .or. min_ref > min_loc ) then
          write(*,*) 'xxx ERROR: REQUESTED DOMAIN IS TOO MUCH BROAD'
          write(*,*) 'xxx -- LATITUDINAL direction over the limit'
          write(*,*) 'xxx -- reference max: ', max_ref
          write(*,*) 'xxx -- reference min: ', min_ref
          write(*,*) 'xxx --     local max: ', max_loc
          write(*,*) 'xxx --     local min: ', min_loc
          call PRC_MPIstop
       endif
    endif

    if ( do_zdirec ) then
       max_ref = maxval( lev_org(:,:,:) )
       !min_ref = minval( lev_org(:,:,:) )
       max_loc = maxval( lev_loc(KS-1:KE,:,:) ) ! HALO + 1
       !min_loc = minval( lev_loc(3:KA,:,:) ) ! HALO + 1

       if ( max_ref < max_loc ) then
       !if ( max_ref < max_loc .or. min_ref > min_loc ) then
          write(*,*) 'xxx ERROR: REQUESTED DOMAIN IS TOO MUCH BROAD'
          write(*,*) 'xxx -- VERTICAL direction over the limit'
          write(*,*) 'xxx -- reference max: ', max_ref
          !write(*,*) 'xxx -- reference min: ', min_ref
          write(*,*) 'xxx --     local max: ', max_loc
          !write(*,*) 'xxx --     local min: ', min_loc
          call PRC_MPIstop
       endif
    endif

    return
  end subroutine check_domain_compatibility

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
      basename    ) ! (in)
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none
    real(RP), intent(out) :: u_latlon(KA,IA,JA)
    real(RP), intent(out) :: v_latlon(KA,IA,JA)
    real(RP), intent(in ) :: u_on_map(KA,IA,JA)
    real(RP), intent(in ) :: v_on_map(KA,IA,JA)
    real(RP), intent(in ) :: xlon(IA,JA)
    real(RP), intent(in ) :: xlat(IA,JA)

    character(LEN=*), intent( in) :: basename

    real(RP) :: truelat1, truelat2
    real(RP) :: stand_lon
    real(RP) :: diff
    real(RP) :: alpha
    real(RP) :: sine(IA,JA)
    real(RP) :: cose(IA,JA)
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

    do j = 1, JA
    do i = 1, IA
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

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       u_latlon(k,i,j) = v_on_map(k,i,j)*sine(i,j) + u_on_map(k,i,j)*cose(i,j)
       v_latlon(k,i,j) = v_on_map(k,i,j)*cose(i,j) - u_on_map(k,i,j)*sine(i,j)
    enddo
    enddo
    enddo

    return
  end subroutine wrf_arwpost_calc_uvmet

  !-----------------------------------------------------------------------------
  ! Haversine Formula (from R.W. Sinnott, "Virtues of the Haversine",
  ! Sky and Telescope, vol. 68, no. 2, 1984, p. 159):
  function haversine( &
      la0,       &
      lo0,       &
      la,        &
      lo )       &
      result( d )
    implicit none
    real(RP), intent(in) :: la0, lo0, la, lo   ! la,la0: Lat, lo,lo0: Lon; [rad]
    real(RP) :: d, dlon, dlat, work1, work2
    !---------------------------------------------------------------------------

    ! output unit : [m]
    dlon = lo0 - lo
    dlat = la0 - la
    work1 = (sin(dlat/2.0_RP))**2.0_RP + &
            cos(la0) * cos(la) * (sin(dlon/2.0_RP))**2.0_RP
    work2 = 2.0_RP * asin(min( 1.0_RP, sqrt(work1) ))
    d = r_in_m * work2

  end function haversine

end module mod_realinput
!-------------------------------------------------------------------------------
