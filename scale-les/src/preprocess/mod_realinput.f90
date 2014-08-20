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
     lon   => REAL_LON,   &
     lat   => REAL_LAT,   &
     lon_u => REAL_LONX,  &
     lat_v => REAL_LATY,  &
     cz    => REAL_CZ,    &
     fz    => REAL_FZ
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
  use scale_const, only:    &
     pi     => CONST_PI,     &
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
     I_SW   => CONST_I_SW,   &
     I_LW   => CONST_I_LW,   &
     CONST_UNDEF
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
  private :: InputSurfaceSCALE
  private :: InputSurfaceWRF
  private :: LatLonZ_Interpolation_Linear
  private :: latlonz_interporation_fact
  private :: haversine

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: iSCALE  = 1
  integer, private, parameter :: iWRFARW = 2
  integer, private, parameter :: iJMAMSM = 3
  integer, private, parameter :: iNICAM  = 4

  real(RP), private, parameter :: large_number_one   = 9.999E+15_RP
  real(RP), private, parameter :: large_number_two   = 8.888E+15_RP
  real(RP), private, parameter :: large_number_three = 7.777E+15_RP
  integer, parameter :: itp_nh = 3
  integer, parameter :: itp_nv = 2

  integer, private :: interp_search_divnum = 10
  logical, private :: wrfout = .false.  ! file type switch (wrfout or wrfrst)

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

    intrinsic shape
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
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
       call ExternalFileGetShape( dims(:), timelen, mdlid,        &
                                   BASENAME_ORG, myrank, single=.true. )
       if ( wrf_file_type_in ) then
          wrfout = .true.
          if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF History Output'
       else
          wrfout = .false.
          if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF Restart'
       endif

    !case('JMA-MSM')

    !case('NCEP-FNL')

    !case('NICAM')

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
      dens,           &
      momz,           &
      momx,           &
      momy,           &
      rhot,           &
      qtrc,           &
      basename_org,   &
      dims,           &
      mdlid,          &
      mptype_parent,  &    
      timelen,        &
      serial          )
    implicit none

    real(RP),         intent(out) :: dens(:,:,:,:)
    real(RP),         intent(out) :: momz(:,:,:,:)
    real(RP),         intent(out) :: momx(:,:,:,:)
    real(RP),         intent(out) :: momy(:,:,:,:)
    real(RP),         intent(out) :: rhot(:,:,:,:)
    real(RP),         intent(out) :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in)  :: basename_org
    integer,          intent(in)  :: dims(:)
    integer,          intent(in)  :: mdlid            ! model type id
    character(LEN=*), intent(in)  :: mptype_parent    ! microphysics type of the parent model (single or double)
    integer,          intent(in)  :: timelen          ! time steps in one file
    logical,          intent(in)  :: serial           ! read by a serial process

    real(RP) :: phy_rd_init_value = 0.0_RP

    integer :: k, i, j
    character(len=H_SHORT)  :: mptype_run

    intrinsic shape
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
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
                            1,                   &  !suffix of start time: ts
                            timelen              )  !suffix of end time: te

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

    integer :: k, i, j, n
    integer :: ts, te
    integer :: timetarg
    !---------------------------------------------------------------------------

    ts = 1
    te = numsteps

    allocate( buffer(KA,IA,JA,te-ts+1) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Boundary]'

    call FILEIO_write( DENS(:,:,:,ts:te),                            &
                       basename, title,                              &
                       'DENS', 'Reference Density', 'kg/m3', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt           )

    do n = ts, te
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMZ(k,i,j,n) / ( DENS(k+1,i,j,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer,                                        &
                       basename, title,                               &
                       'VELZ', 'Reference Velocity w', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt            )

    do n = ts, te
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMX(k,i,j,n) / ( DENS(k,i+1,j,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer,                                        &
                       basename, title,                               &
                       'VELX', 'Reference Velocity u', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt            )

    do n = ts, te
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMY(k,i,j,n) / ( DENS(k,i,j+1,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer,                                        &
                       basename, title,                               &
                       'VELY', 'Reference Velocity y', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt            )

    do n = ts, te
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = RHOT(k,i,j,n) / DENS(k,i,j,n)
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer,                                        &
                       basename, title,                               &
                       'POTT', 'Reference PT', 'K', 'ZXYT',           &
                       atmos_boundary_out_dtype, update_dt            )

    do n = ts, te
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = QTRC(k,i,j,n,I_QV)
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer,                                         &
                       basename, title,                                &
                       'QV', 'Reference water vapor', 'kg/kg', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt             )

    deallocate( buffer )

    return
  end subroutine ParentAtomBoundary

  !-----------------------------------------------------------------------------
  !> Land&Ocean Data Read/Write
  !-----------------------------------------------------------------------------
  subroutine ParentSurfaceInput( &
      dens,         &
      momz,         &
      momx,         &
      momy,         &
      rhot,         &
      qtrc,         &
      basename_org, &
      dims,         &
      timelen,      &
      mdlid,        &
      serial,       &
      no_add_input  )
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
    integer,          intent(in) :: timelen  ! time steps in one file
    integer,          intent(in) :: mdlid    ! model type id
    logical,          intent(in) :: serial   ! read by a serial process
    logical,          intent(in) :: no_add_input

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
    real(RP) :: albw (IA,JA,2)
    real(RP) :: albg (IA,JA,2)
    real(RP) :: z0w  (IA,JA)
    real(RP) :: skint(IA,JA)
    real(RP) :: skinw(IA,JA)
    real(RP) :: snowq(IA,JA)
    real(RP) :: snowt(IA,JA)

    real(RP) :: tc_urb(IA,JA)
    real(RP) :: qc_urb(IA,JA)
    real(RP) :: uc_urb(IA,JA)

    integer :: i, j

    intrinsic shape
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Input-Surface]'

    if( LKMAX < 4 )then
       write(*,*) 'xxx LKMAX less than 4: ', LKMAX
       write(*,*) 'xxx in Real Case, LKMAX should be set more than 4'
       call PRC_MPIstop
    endif

    ! Read Data
    if( mdlid == iSCALE ) then ! TYPE: SCALE-LES
       call InputSurfaceSCALE( tg   (:,:,:), &
                               strg (:,:,:), &
                               roff (:,:),   &
                               qvef (:,:),   &
                               tw   (:,:),   &
                               lst  (:,:),   &
                               ust  (:,:),   &
                               sst  (:,:),   &
                               albw (:,:,:), &
                               albg (:,:,:), &
                               z0w  (:,:),   &
                               skint(:,:),   &
                               skinw(:,:),   &
                               snowq(:,:),   &
                               snowt(:,:),   &
                               basename_org  )
    elseif( mdlid == iWRFARW ) then ! TYPE: WRF-ARW
       call InputSurfaceWRF( tg   (:,:,:), &
                             strg (:,:,:), &
                             roff (:,:),   &
                             qvef (:,:),   &
                             tw   (:,:),   &
                             lst  (:,:),   &
                             ust  (:,:),   &
                             sst  (:,:),   &
                             albw (:,:,:), &
                             albg (:,:,:), &
                             z0w  (:,:),   &
                             skint(:,:),   &
                             skinw(:,:),   &
                             snowq(:,:),   &
                             snowt(:,:),   &
                             basename_org, &
                             dims(:),      &
                             mdlid,        &
                             serial,       &
                             no_add_input  )
    endif

    call THERMODYN_temp_pres( temp(:,:,:),     & ! [OUT]
                              pres(:,:,:),     & ! [OUT]
                              dens(:,:,:,1),   & ! [IN]
                              rhot(:,:,:,1),   & ! [IN]
                              qtrc(:,:,:,:,1)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
       tc_urb(i,j) = temp(KS,i,j)
       qc_urb(i,j) = qtrc(KS,i,j,1,I_QV)
       uc_urb(i,j) = sqrt( ( momx(KS,i,j,1) / ( dens(KS,i+1,  j,1) + dens(KS,i,j,1) ) * 2.0_RP )**2.0_RP &
                         + ( momy(KS,i,j,1) / ( dens(KS,  i,j+1,1) + dens(KS,i,j,1) ) * 2.0_RP )**2.0_RP )
    enddo
    enddo

    ! Write out: Ocean
    call OCEAN_vars_external_in( tw, sst, albw, z0w )

    ! Write out: Land
    call LAND_vars_external_in( tg, strg, lst, albg )

    ! Write out: Urban
    call URBAN_vars_external_in( ust, tc_urb, qc_urb, uc_urb )

    ! Input to PHY_SF Container
    call ATMOS_PHY_SF_vars_external_in( skint, albw(:,:,I_LW), albw(:,:,I_SW), z0w )

    return
  end subroutine ParentSurfaceInput

  !-----------------------------------------------------------------------------
  !> Individual Procedures
  !-----------------------------------------------------------------------------
  subroutine InputAtomSCALE( &
      dens,          & ! (out)
      momz,          & ! (out)
      momx,          & ! (out)
      momy,          & ! (out)
      rhot,          & ! (out)
      qtrc,          & ! (out)
      basename_org,  & ! (in)
      start_step,    & ! (in)
      end_step       ) ! (in)
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
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    implicit none

    real(RP),         intent(out) :: dens(:,:,:,:)
    real(RP),         intent(out) :: momz(:,:,:,:)
    real(RP),         intent(out) :: momx(:,:,:,:)
    real(RP),         intent(out) :: momy(:,:,:,:)
    real(RP),         intent(out) :: rhot(:,:,:,:)
    real(RP),         intent(out) :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in)  :: basename_org
    integer,          intent(in)  :: start_step
    integer,          intent(in)  :: end_step

    ! work
    real(RP), allocatable :: read2D(:,:)
    real(RP), allocatable :: read3D(:,:,:)

    real(RP), allocatable :: LON_ORG (:,:,:)
    real(RP), allocatable :: LAT_ORG (:,:,:)
    real(RP), allocatable :: LONU_ORG(:,:,:)
    real(RP), allocatable :: LATU_ORG(:,:,:)
    real(RP), allocatable :: LONV_ORG(:,:,:)
    real(RP), allocatable :: LATV_ORG(:,:,:)
    real(RP), allocatable :: GEOH_ORG(:,:,:,:)
    real(RP), allocatable :: GEOF_ORG(:,:,:,:)

    real(RP), allocatable :: DENS_ORG(:,:,:,:)
    real(RP), allocatable :: MOMZ_ORG(:,:,:,:)
    real(RP), allocatable :: MOMX_ORG(:,:,:,:)
    real(RP), allocatable :: MOMY_ORG(:,:,:,:)
    real(RP), allocatable :: RHOT_ORG(:,:,:,:)
    real(RP), allocatable :: POTT_ORG(:,:,:,:)
    real(RP), allocatable :: QTRC_ORG(:,:,:,:,:)
    real(RP), allocatable :: PSFC_ORG(:,:,:)
    real(RP), allocatable :: TSFC_ORG(:,:,:)

    real(RP), allocatable :: W_ORG(:,:,:,:)
    real(RP), allocatable :: U_ORG(:,:,:,:)
    real(RP), allocatable :: V_ORG(:,:,:,:)

    real(RP) :: W   (KA,IA,JA,start_step:end_step)
    real(RP) :: U   (KA,IA,JA,start_step:end_step)
    real(RP) :: V   (KA,IA,JA,start_step:end_step)
    real(RP) :: POTT(KA,IA,JA,start_step:end_step)

    real(RP) :: PRES(KA,IA,JA,start_step:end_step)
    real(RP) :: TEMP(KA,IA,JA,start_step:end_step)
    real(RP) :: QV  (KA,IA,JA,start_step:end_step)
    real(RP) :: QC  (KA,IA,JA,start_step:end_step)

    real(RP) :: temp_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: pott_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: pres_sfc(1,IA,JA,start_step:end_step)
    real(RP) :: qv_sfc  (1,IA,JA,start_step:end_step)
    real(RP) :: qc_sfc  (1,IA,JA,start_step:end_step)

    integer :: dims(6) ! 1-3: cell center [z,x,y], 4-6: staggered [z,x,y]
    integer :: KALL, IALL, JALL
    integer :: rank

    integer :: k, i, j, n, iq
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye

    !logical :: single_ = .false.

    intrinsic shape
    !---------------------------------------------------------------------------

    KALL = PARENT_KMAX
    IALL = PARENT_IMAX * NEST_TILE_NUM_X
    JALL = PARENT_JMAX * NEST_TILE_NUM_Y

    allocate( read2D( PARENT_IMAX, PARENT_JMAX              ) )
    allocate( read3D( PARENT_IMAX, PARENT_JMAX, PARENT_KMAX ) )

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
    allocate( pott_org( KALL, IALL, JALL, start_step:end_step )    )
    allocate( qtrc_org( KALL, IALL, JALL, start_step:end_step, QA) )
    allocate( psfc_org(       IALL, JALL, start_step:end_step )    )
    allocate( tsfc_org(       IALL, JALL, start_step:end_step )    )

    allocate( w_org   ( KALL, IALL, JALL, start_step:end_step )    )
    allocate( u_org   ( KALL, IALL, JALL, start_step:end_step )    )
    allocate( v_org   ( KALL, IALL, JALL, start_step:end_step )    )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputSCALE]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX * (xloc-1) + 1
       xe = PARENT_IMAX * xloc
       ys = PARENT_JMAX * (yloc-1) + 1
       ye = PARENT_JMAX * yloc

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

       do n = start_step, end_step
         call FileRead( read2D(:,:), BASENAME_ORG, "SFC_PRES", n, rank )
         psfc_org(xs:xe,ys:ye,n) = read2D(:,:)
       end do

       do n = start_step, end_step
         call FileRead( read2D(:,:), BASENAME_ORG, "SFC_TEMP", n, rank )
         tsfc_org(xs:xe,ys:ye,n) = read2D(:,:)
       end do

    end do

    ! make dimension
    dims(1) = KALL ! z center
    dims(2) = IALL ! x center
    dims(3) = JALL ! y center
    dims(4) = KALL ! z staggered
    dims(5) = IALL ! x staggered
    dims(6) = JALL ! y staggered

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

    ! convert from momentum to velocity
    do n = start_step, end_step
    do j = 1, JALL
    do i = 1, IALL
    do k = 1, KALL-1
       w_org(k,i,j,n) = 2.0_RP * momz_org(k,i,j,n) / ( dens_org(k+1,i,j,n) + dens_org(k,i,j,n) )
    end do
    end do
    end do
    end do
    do n = start_step, end_step
    do j = 1, JALL
    do i = 1, IALL
       w_org(KALL,i,j,n) = 0.0_RP
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JALL
    do i = 1, IALL-1
    do k = 1, KALL
       u_org(k,i,j,n) = 2.0_RP * momx_org(k,i,j,n) / ( dens_org(k,i+1,j,n) + dens_org(k,i,j,n) )
    end do
    end do
    end do
    end do
    do n = start_step, end_step
    do j = 1, JALL
    do k = 1, KALL
       u_org(k,IALL,j,n) = momx_org(k,IALL,j,n) / dens_org(k,IALL,j,n)
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JALL-1
    do i = 1, IALL
    do k = 1, KALL
       v_org(k,i,j,n) = 2.0_RP * momy_org(k,i,j,n) / ( dens_org(k,i,j+1,n) + dens_org(k,i,j,n) )
    end do
    end do
    end do
    end do
    do n = start_step, end_step
    do i = 1, IALL
    do k = 1, KALL
       v_org(k,i,JALL,n) = momy_org(k,i,JALL,n) / dens_org(k,i,JALL,n)
    end do
    end do
    end do

    do n = start_step, end_step
    do j = 1, JALL
    do i = 1, IALL
    do k = 1, KALL
       pott_org(k,i,j,n) = rhot_org(k,i,j,n) / dens_org(k,i,j,n)
       dens_org(k,i,j,n) = log( dens_org(k,i,j,n) )
    end do
    end do
    end do
    end do

    ! make initial condition by interpolation
    call LatLonZ_Interpolation_Linear( dens    (:,:,:,:),   &
                                       w       (:,:,:,:),   &
                                       u       (:,:,:,:),   &
                                       v       (:,:,:,:),   &
                                       pott    (:,:,:,:),   &
                                       qtrc    (:,:,:,:,:), &
                                       pres_sfc(1,:,:,:),   &
                                       pott_sfc(1,:,:,:),   &
                                       dens_org(:,:,:,:),   &
                                       w_org   (:,:,:,:),   &
                                       u_org   (:,:,:,:),   &
                                       v_org   (:,:,:,:),   &
                                       pott_org(:,:,:,:),   &
                                       qtrc_org(:,:,:,:,:), &
                                       psfc_org(:,:,:),     &
                                       tsfc_org(:,:,:),     &
                                       lat_org (:,:,:),     &
                                       lon_org (:,:,:),     &
                                       latu_org(:,:,:),     &
                                       lonu_org(:,:,:),     &
                                       latv_org(:,:,:),     &
                                       lonv_org(:,:,:),     &
                                       geof_org(:,:,:,:),   &
                                       geoh_org(:,:,:,:),   &
                                       dims    (:),         &
                                       start_step,          &
                                       end_step             )

    do n = start_step, end_step !--- time loop

      if ( I_QV > 0 ) then
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            qv(k,i,j,n) = qtrc(k,i,j,n,I_QV)
         end do
         end do
         end do
      else
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            qv(k,i,j,n) = 0.0_RP
         end do
         end do
         end do
      end if

      if ( I_QC > 0 ) then
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            qc(k,i,j,n) = qtrc(k,i,j,n,I_QC)
         end do
         end do
         end do
      else
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            qc(k,i,j,n) = 0.0_RP
         end do
         end do
         end do
      end if

      ! make density & pressure profile in moist condition
      qv_sfc(1,:,:,n) = qv(KS,:,:,n)
      qc_sfc(1,:,:,n) = qc(KS,:,:,n)
      call HYDROSTATIC_buildrho( dens    (:,:,:,n), & ! [OUT]
                                 temp    (:,:,:,n), & ! [OUT]
                                 pres    (:,:,:,n), & ! [OUT]
                                 pott    (:,:,:,n), & ! [IN]
                                 qv      (:,:,:,n), & ! [IN]
                                 qc      (:,:,:,n), & ! [IN]
                                 temp_sfc(:,:,:,n), & ! [OUT]
                                 pres_sfc(:,:,:,n), & ! [IN]
                                 pott_sfc(:,:,:,n), & ! [IN]
                                 qv_sfc  (:,:,:,n), & ! [IN]
                                 qc_sfc  (:,:,:,n)  ) ! [IN]

      call COMM_vars8( dens(:,:,:,n), 1 )
      call COMM_wait ( dens(:,:,:,n), 1 )

      do j = JS, JE
      do i = IS, IE
      do k = KS, KE
         momz(k,i,j,n) = 0.5_RP * w(k,i,j,n) * ( dens(k+1,  i,  j,n) + dens(k,i,j,n) )
         momx(k,i,j,n) = 0.5_RP * u(k,i,j,n) * ( dens(k,  i+1,  j,n) + dens(k,i,j,n) )
         momy(k,i,j,n) = 0.5_RP * v(k,i,j,n) * ( dens(k,    i,j+1,n) + dens(k,i,j,n) )
         rhot(k,i,j,n) = pott(k,i,j,n) * dens(k,i,j,n)
      enddo
      enddo
      enddo

    enddo !--- time loop

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
    deallocate( pott_org )
    deallocate( qtrc_org )
    deallocate( psfc_org )
    deallocate( tsfc_org )

    deallocate( w_org    )
    deallocate( u_org    )
    deallocate( v_org    )

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
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    implicit none

    real(RP),         intent(out)  :: dens(:,:,:,:)
    real(RP),         intent(out)  :: momz(:,:,:,:)
    real(RP),         intent(out)  :: momx(:,:,:,:)
    real(RP),         intent(out)  :: momy(:,:,:,:)
    real(RP),         intent(out)  :: rhot(:,:,:,:)
    real(RP),         intent(out)  :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in)   :: basename
    character(LEN=*), intent(in)   :: mptype_run
    character(LEN=*), intent(in)   :: mptype_parent
    integer,          intent(in)   :: mdlid
    integer,          intent(in)   :: dims(:)
    integer,          intent(in)   :: ts      ! suffix of start time
    integer,          intent(in)   :: te      ! suffix of end time
    logical,          intent(in)   :: serial  ! read by a serial process

    real(RP), allocatable :: DENS_ORG(:,:,:,:)
    real(RP), allocatable :: P_ORG(:,:,:,:)
    real(RP), allocatable :: PBASE_ORG(:,:,:,:)
    real(RP), allocatable :: W_ORG(:,:,:,:)
    real(RP), allocatable :: U_ORG(:,:,:,:)
    real(RP), allocatable :: V_ORG(:,:,:,:)
    real(RP), allocatable :: POTT_ORG(:,:,:,:)
    real(RP), allocatable :: QTRC_ORG(:,:,:,:,:)
    real(RP), allocatable :: LAT_ORG(:,:,:)
    real(RP), allocatable :: LON_ORG(:,:,:)
    real(RP), allocatable :: LATU_ORG(:,:,:)
    real(RP), allocatable :: LONU_ORG(:,:,:)
    real(RP), allocatable :: LATV_ORG(:,:,:)
    real(RP), allocatable :: LONV_ORG(:,:,:)
    real(RP), allocatable :: GEOF_ORG(:,:,:,:)
    real(RP), allocatable :: GEOH_ORG(:,:,:,:)
    real(RP), allocatable :: PSFC_ORG(:,:,:)
    real(RP), allocatable :: TH2M_ORG(:,:,:)

    real(RP), allocatable :: AQ_CV(:) !< CV for each hydrometeors [J/kg/K]

    real(SP), allocatable :: dummy_3D(:,:,:)     ! for WRF restart file
    real(SP), allocatable :: dummy_4D(:,:,:,:)   ! for WRF restart file
    real(SP), allocatable :: dummy_5D(:,:,:,:,:) ! for WRF restart file

    real(RP) :: W   (KA,IA,JA,te)
    real(RP) :: U   (KA,IA,JA,te)
    real(RP) :: V   (KA,IA,JA,te)
    real(RP) :: POTT(KA,IA,JA,te)

    real(RP) :: PRES(KA,IA,JA,te)
    real(RP) :: TEMP(KA,IA,JA,te)
    real(RP) :: QV  (KA,IA,JA,te)
    real(RP) :: QC  (KA,IA,JA,te)

    real(RP) :: temp_sfc(1,IA,JA,te)
    real(RP) :: pott_sfc(1,IA,JA,te)
    real(RP) :: pres_sfc(1,IA,JA,te)
    real(RP) :: qv_sfc  (1,IA,JA,te)
    real(RP) :: qc_sfc  (1,IA,JA,te)

    real(RP) :: qdry
    real(RP) :: CVtot
    real(RP) :: Rtot
    real(RP) :: RovCP
    real(RP) :: tem
    real(RP) :: t0 = 300.0_RP  ! Defined in WRF

    integer :: n, k, i, j, iq, iqw, iq_all
    character(len=H_MID)  :: varname_wrf
    logical :: do_read

    intrinsic shape
    !---------------------------------------------------------------------------

    allocate( AQ_CV(QQA) )
    AQ_CV(I_QV) = CVvap
    if ( QWS /= 0 ) then
       do n = QWS, QWE
          AQ_CV(n) = CL
       enddo
    endif
    if ( QIS /= 0 ) then
       do n = QIS, QIE
          AQ_CV(n) = CI
       enddo
    endif

    allocate( p_org(dims(1),dims(2),dims(3),te) )
    allocate( pbase_org(dims(1),dims(2),dims(3),te) )
    allocate( dens_org(dims(1),dims(2),dims(3),te) )
    allocate( w_org(dims(4),dims(2),dims(3),te) )
    allocate( u_org(dims(1),dims(5),dims(3),te) )
    allocate( v_org(dims(1),dims(2),dims(6),te) )
    allocate( pott_org(dims(1),dims(2),dims(3),te) )
    allocate( qtrc_org(dims(1),dims(2),dims(3),te,QA) )
    allocate( psfc_org(dims(2),dims(3),te) )
    allocate( th2m_org(dims(2),dims(3),te) )
    allocate( lat_org(dims(2),dims(3),te) )
    allocate( lon_org(dims(2),dims(3),te) )
    allocate( latu_org(dims(5),dims(3),te) )
    allocate( lonu_org(dims(5),dims(3),te) )
    allocate( latv_org(dims(2),dims(6),te) )
    allocate( lonv_org(dims(2),dims(6),te) )
    allocate( geoh_org(dims(1),dims(2),dims(3),te) )
    allocate( geof_org(dims(4),dims(2),dims(3),te) )

    if( serial ) then
       if( myrank == PRC_master ) then
          do_read = .true.
       else
          do_read = .false.
       endif
    else
       do_read = .true.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputWRF]'

    allocate( dummy_4D(dims(1),dims(2),dims(3),te) )
    if( do_read ) then
       call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                      BASENAME, "P",  ts, te, myrank, mdlid, single=.true. )
       p_org(:,:,:,:) = dummy_4D(:,:,:,:)
       call ExternalFileRead( dummy_4D(:,:,:,:),                         &
                      BASENAME, "PB", ts, te, myrank, mdlid, single=.true. )
       pbase_org(:,:,:,:) = dummy_4D(:,:,:,:)
    endif

    if ( wrfout ) then
       varname_wrf = "T"
    else
       varname_wrf = "T_1"
    endif
    if( do_read ) then
       call ExternalFileRead( dummy_4D(:,:,:,:),                          &
                      BASENAME, varname_wrf,  ts, te, myrank, mdlid, single=.true. )
       pott_org(:,:,:,:) = dummy_4D(:,:,:,:)
    endif
    deallocate( dummy_4D )

    allocate( dummy_4D(dims(4),dims(2),dims(3),te) )
    if ( wrfout ) then
       varname_wrf = "W"
    else
       varname_wrf = "W_1"
    endif
    if( do_read ) then
       call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                      BASENAME, varname_wrf,  ts, te, myrank, mdlid, single=.true., zstag=.true. )
       w_org(:,:,:,:) = dummy_4D(:,:,:,:)
    endif
    deallocate( dummy_4D )
    allocate( dummy_4D(dims(1),dims(5),dims(3),te) )
    if ( wrfout ) then
       varname_wrf = "U"
    else
       varname_wrf = "U_1"
    endif
    if( do_read ) then
       call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                      BASENAME, varname_wrf,  ts, te, myrank, mdlid, single=.true., xstag=.true. )
       u_org(:,:,:,:) = dummy_4D(:,:,:,:)
    endif
    deallocate( dummy_4D )
    allocate( dummy_4D(dims(1),dims(2),dims(6),te) )
    if ( wrfout ) then
       varname_wrf = "V"
    else
       varname_wrf = "V_1"
    endif
    if( do_read ) then
       call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                      BASENAME, varname_wrf,  ts, te, myrank, mdlid, single=.true., ystag=.true. )
       v_org(:,:,:,:) = dummy_4D(:,:,:,:)
    endif
    deallocate( dummy_4D )

    allocate( dummy_5D(dims(1),dims(2),dims(3),te,QA) )
    dummy_5D(:,:,:,:,:) = 0.0_SP
    if( do_read ) then
       call ExternalFileRead( dummy_5D(:,:,:,:,I_QV),                     &
                      BASENAME, "QVAPOR", ts, te, myrank, mdlid, single=.true. )
    endif

    if( trim(mptype_run)/='dry' )then
       if( do_read ) then
          call ExternalFileRead( dummy_5D(:,:,:,:,I_QC),                     &
                         BASENAME, "QCLOUD", ts, te, myrank, mdlid, single=.true. )
          call ExternalFileRead( dummy_5D(:,:,:,:,I_QR),                     &
                         BASENAME, "QRAIN",  ts, te, myrank, mdlid, single=.true. )
          call ExternalFileRead( dummy_5D(:,:,:,:,I_QI),                     &
                         BASENAME, "QICE",   ts, te, myrank, mdlid, single=.true. )
          call ExternalFileRead( dummy_5D(:,:,:,:,I_QS),                     &
                         BASENAME, "QSNOW",  ts, te, myrank, mdlid, single=.true. )
          call ExternalFileRead( dummy_5D(:,:,:,:,I_QG),                     &
                         BASENAME, "QGRAUP", ts, te, myrank, mdlid, single=.true. )
       endif
    endif

    if( trim(mptype_run)=='dry' )then
       iq_all = I_QV
    elseif( trim(mptype_run)=='single' )then
       iq_all = I_QG
    elseif( trim(mptype_run)=='double' )then
       iq_all = I_NG

       if( trim(mptype_parent)=='double' )then
          if( do_read ) then
             call ExternalFileRead( dummy_5D(:,:,:,:,I_NC),                     &
                            BASENAME, "NC", ts, te, myrank, mdlid, single=.true. )
             call ExternalFileRead( dummy_5D(:,:,:,:,I_NR),                     &
                            BASENAME, "NR", ts, te, myrank, mdlid, single=.true. )
             call ExternalFileRead( dummy_5D(:,:,:,:,I_NI),                     &
                            BASENAME, "NI", ts, te, myrank, mdlid, single=.true. )
             call ExternalFileRead( dummy_5D(:,:,:,:,I_NS),                     &
                            BASENAME, "NS", ts, te, myrank, mdlid, single=.true. )
             call ExternalFileRead( dummy_5D(:,:,:,:,I_NG),                     &
                            BASENAME, "NG", ts, te, myrank, mdlid, single=.true. )
          endif
       endif
    endif

    qtrc_org(:,:,:,:,:) = dummy_5D(:,:,:,:,:)
    deallocate( dummy_5D )

    do iq = I_QV, iq_all
    do n = ts, te
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       qtrc_org(k,i,j,n,iq) = max( qtrc_org(k,i,j,n,iq), 0.0_RP )
    enddo
    enddo
    enddo
    enddo
    enddo

    if( trim(mptype_run)=='double' .and. trim(mptype_parent)=='single' )then
       if( IO_L ) write(IO_FID_LOG,*) '--- Diagnose Number Concentration from Mixing Ratio'
       call diagnose_number_concentration( qtrc_org(:,:,:,:,:) ) ! [inout]
    endif


    allocate( dummy_3D(dims(2),dims(3),te) )
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "XLAT",    ts, te, myrank, mdlid, single=.true. )
       lat_org(:,:,:) = dummy_3D(:,:,:) * d2r
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "XLONG",   ts, te, myrank, mdlid, single=.true. )
       lon_org(:,:,:) = dummy_3D(:,:,:) * d2r
    endif
    deallocate( dummy_3D )
    allocate( dummy_3D(dims(5),dims(3),te) )
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "XLAT_U",  ts, te, myrank, mdlid, single=.true., xstag=.true. )
       latu_org(:,:,:) = dummy_3D(:,:,:) * d2r
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "XLONG_U", ts, te, myrank, mdlid, single=.true., xstag=.true. )
       lonu_org(:,:,:) = dummy_3D(:,:,:) * d2r
    endif
    deallocate( dummy_3D )
    allocate( dummy_3D(dims(2),dims(6),te) )
    if( do_read ) then
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLAT_V",  ts, te, myrank, mdlid, single=.true., ystag=.true. )
    latv_org(:,:,:) = dummy_3D(:,:,:) * d2r
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLONG_V", ts, te, myrank, mdlid, single=.true., ystag=.true. )
    lonv_org(:,:,:) = dummy_3D(:,:,:) * d2r
    endif
    deallocate( dummy_3D )

    allocate( dummy_4D(dims(4),dims(2),dims(3),te) )
    if( do_read ) then
       call ExternalFileRead( dummy_4D(:,:,:,:),                           &
                      BASENAME, "PHB", ts, te, myrank, mdlid, single=.true., zstag=.true. )
       geof_org(:,:,:,:) = dummy_4D(:,:,:,:)
       call ExternalFileRead( dummy_4D(:,:,:,:),                           &
                      BASENAME, "PH",  ts, te, myrank, mdlid, single=.true., zstag=.true. )
       geof_org(:,:,:,:) = geof_org(:,:,:,:) + dummy_4D(:,:,:,:)
    endif
    deallocate( dummy_4D )

    allocate( dummy_3D(dims(2),dims(3),te) )
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                           &
                      BASENAME, "PSFC", ts, te, myrank, mdlid, single=.true. )
       psfc_org(:,:,:) = dummy_3D(:,:,:)
       call ExternalFileRead( dummy_3D(:,:,:),                           &
                      BASENAME, "TH2",  ts, te, myrank, mdlid, single=.true. )
       th2m_org(:,:,:) = dummy_3D(:,:,:)
    endif
    deallocate( dummy_3D )

    if( serial ) then
       call COMM_bcast( p_org, dims(1),dims(2),dims(3),te )
       call COMM_bcast( pbase_org, dims(1),dims(2),dims(3),te )
       call COMM_bcast( pott_org, dims(1),dims(2),dims(3),te )
       call COMM_bcast( w_org, dims(4),dims(2),dims(3),te )
       call COMM_bcast( u_org, dims(1),dims(5),dims(3),te )
       call COMM_bcast( v_org, dims(1),dims(2),dims(6),te )

       do iq = I_QV, iq_all
          call COMM_bcast( qtrc_org(:,:,:,:,iq), dims(1),dims(2),dims(3),te )
       enddo

       call COMM_bcast( lat_org(:,:,:), dims(2),dims(3),te )
       call COMM_bcast( lon_org(:,:,:), dims(2),dims(3),te )
       call COMM_bcast( latu_org(:,:,:), dims(5),dims(3),te )
       call COMM_bcast( lonu_org(:,:,:), dims(5),dims(3),te )
       call COMM_bcast( latv_org(:,:,:), dims(2),dims(6),te )
       call COMM_bcast( lonv_org(:,:,:), dims(2),dims(6),te )

       call COMM_bcast( geof_org(:,:,:,:), dims(4),dims(2),dims(3),te )
       call COMM_bcast( psfc_org(:,:,:), dims(2),dims(3),te )
       call COMM_bcast( th2m_org(:,:,:), dims(2),dims(3),te )
    endif

    do n = 1, te
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       p_org(k,i,j,n) = pbase_org(k,i,j,n) + p_org(k,i,j,n)
       pott_org(k,i,j,n) = pott_org(k,i,j,n) + t0
    end do
    end do
    end do
    end do

    ! convert to geopotential height to use as real height in WRF
    do n = 1, te
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       geof_org(k,i,j,n) = geof_org(k,i,j,n) / grav
    end do
    end do
    end do
    end do
    ! make half level of geopotential height from face level
    do n = 1, te
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       geoh_org(k,i,j,n) = (geof_org(k,i,j,n) + geof_org(k+1,i,j,n))*0.5_RP
    end do
    end do
    end do
    end do

    ! calc dens
    do n = ts, te
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       qdry  = 1.0_RP
       CVtot = 0.0_RP
       do iqw = QQS, QQE
          qdry  = qdry  - qtrc_org(k,i,j,n,iqw)
          CVtot = CVtot + qtrc_org(k,i,j,n,iqw) * AQ_CV(iqw)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * qtrc_org(k,i,j,n,I_QV)
       RovCP = Rtot / ( CVtot + Rtot )
       tem   = pott_org(k,i,j,n) * ( p_org(k,i,j,n) / PRE00 )**RovCP

       dens_org(k,i,j,n) = p_org(k,i,j,n) / ( Rtot * tem )
    end do
    end do
    end do
    end do

    ! make initial condition by interpolation
    call LatLonZ_Interpolation_Linear( dens,w,u,v,pott,qtrc,pres_sfc(1,:,:,:),pott_sfc(1,:,:,:), &
                                       dens_org,w_org,u_org,v_org,pott_org,qtrc_org,        &
                                       psfc_org,th2m_org,                                   &
                                       lat_org,lon_org,latu_org,lonu_org,latv_org,lonv_org, &
                                       geof_org,geoh_org,dims,ts,te )

    do n = ts, te !--- time loop

      if ( I_QV > 0 ) then
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            qv(k,i,j,n) = qtrc(k,i,j,n,I_QV)
         end do
         end do
         end do
      else
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            qv(k,i,j,n) = 0.0_RP
         end do
         end do
         end do
      end if

      if ( I_QC > 0 ) then
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            qc(k,i,j,n) = qtrc(k,i,j,n,I_QC)
         end do
         end do
         end do
      else
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            qc(k,i,j,n) = 0.0_RP
         end do
         end do
         end do
      end if

      ! make density & pressure profile in moist condition
      qv_sfc(1,:,:,n) = qv(KS,:,:,n)
      qc_sfc(1,:,:,n) = qc(KS,:,:,n)
      call HYDROSTATIC_buildrho( dens    (:,:,:,n), & ! [OUT]
                                 temp    (:,:,:,n), & ! [OUT]
                                 pres    (:,:,:,n), & ! [OUT]
                                 pott    (:,:,:,n), & ! [IN]
                                 qv      (:,:,:,n), & ! [IN]
                                 qc      (:,:,:,n), & ! [IN]
                                 temp_sfc(:,:,:,n), & ! [OUT]
                                 pres_sfc(:,:,:,n), & ! [IN]
                                 pott_sfc(:,:,:,n), & ! [IN]
                                 qv_sfc  (:,:,:,n), & ! [IN]
                                 qc_sfc  (:,:,:,n)  ) ! [IN]

      call COMM_vars8( dens(:,:,:,n), 1 )
      call COMM_wait ( dens(:,:,:,n), 1 )

      do j = JS, JE
      do i = IS, IE
      do k = KS, KE
         momz(k,i,j,n) = 0.5_RP * w(k,i,j,n) * ( dens(k+1,  i,  j,n) + dens(k,i,j,n) )
         momx(k,i,j,n) = 0.5_RP * u(k,i,j,n) * ( dens(k,  i+1,  j,n) + dens(k,i,j,n) )
         momy(k,i,j,n) = 0.5_RP * v(k,i,j,n) * ( dens(k,    i,j+1,n) + dens(k,i,j,n) )
         rhot(k,i,j,n) = pott(k,i,j,n) * dens(k,i,j,n)
      enddo
      enddo
      enddo

    enddo !--- time loop

    deallocate( p_org )
    deallocate( pbase_org )
    deallocate( dens_org )
    deallocate( w_org )
    deallocate( u_org )
    deallocate( v_org )
    deallocate( pott_org )
    deallocate( qtrc_org )

    deallocate( lat_org )
    deallocate( lon_org )
    deallocate( latu_org )
    deallocate( lonu_org )
    deallocate( latv_org )
    deallocate( lonv_org )
    deallocate( geoh_org )
    deallocate( geof_org )
    deallocate( psfc_org )
    deallocate( th2m_org )

    deallocate( AQ_CV )

    return
  end subroutine InputAtomWRF

  !-----------------------------------------------------------------------------
  subroutine InputSurfaceSCALE( &
      tg,          & ! (out)
      strg,        & ! (out)
      roff,        & ! (out)
      qvef,        & ! (out)
      tw,          & ! (out)
      lst,         & ! (out)
      ust,         & ! (out)
      sst,         & ! (out)
      albw,        & ! (out)
      albg,        & ! (out)
      z0w,         & ! (out)
      skint,       & ! (out)
      skinw,       & ! (out)
      snowq,       & ! (out)
      snowt,       & ! (out)
      basename_org ) ! (in)
    use scale_const, only: &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00
    use scale_grid_nest, only: &
       PARENT_KMAX,     &
       PARENT_IMAX,     &
       PARENT_JMAX,     &
       PARENT_LKMAX,    &
       NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y, &
       NEST_TILE_ID
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

    ! work
    real(RP), allocatable :: read1D(:)
    real(RP), allocatable :: read2D(:,:)
    real(RP), allocatable :: read3D(:,:,:)

    real(RP), allocatable :: lon_org(:,:,:)
    real(RP), allocatable :: lat_org(:,:,:)
    real(RP), allocatable :: lz_org (:,:,:,:)

    real(RP), allocatable :: tg_org   (:,:,:)
    real(RP), allocatable :: strg_org (:,:,:)
    real(RP), allocatable :: tw_org   (:,:)
    real(RP), allocatable :: lst_org  (:,:)
    real(RP), allocatable :: ust_org  (:,:)
    real(RP), allocatable :: sst_org  (:,:)
    real(RP), allocatable :: albw_org (:,:,:)
    real(RP), allocatable :: albg_org (:,:,:)
    real(RP), allocatable :: z0w_org  (:,:)
    real(RP), allocatable :: skint_org(:,:)
    real(RP), allocatable :: skinw_org(:,:)
    real(RP), allocatable :: snowq_org(:,:)
    real(RP), allocatable :: snowt_org(:,:)

    real(RP), allocatable :: hfact(:,:,:)
    real(RP), allocatable :: vfact(:,:,:,:,:)

    integer, allocatable :: igrd(:,:,:)
    integer, allocatable :: jgrd(:,:,:)
    integer, allocatable :: kgrd(:,:,:,:,:)

    real(RP) :: lcz_3D(LKMAX,IA,JA)

    integer :: KALL, IALL, JALL
    integer :: rank

    integer :: k, i, j, n
    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye
    !---------------------------------------------------------------------------

    KALL = PARENT_LKMAX
    IALL = PARENT_IMAX * NEST_TILE_NUM_X
    JALL = PARENT_JMAX * NEST_TILE_NUM_Y

    allocate( read1D(                           PARENT_LKMAX ) )
    allocate( read2D( PARENT_IMAX, PARENT_JMAX               ) )
    allocate( read3D( PARENT_IMAX, PARENT_JMAX, PARENT_LKMAX ) )

    allocate( lon_org(       IALL, JALL, 1 ) )
    allocate( lat_org(       IALL, JALL, 1 ) )
    allocate( lz_org ( KALL, IALL, JALL, 1 ) )

    allocate( tg_org   ( KALL, IALL, JALL    ) )
    allocate( strg_org ( KALL, IALL, JALL    ) )
    allocate( tw_org   (       IALL, JALL    ) )
    allocate( lst_org  (       IALL, JALL    ) )
    allocate( ust_org  (       IALL, JALL    ) )
    allocate( sst_org  (       IALL, JALL    ) )
    allocate( albw_org (       IALL, JALL, 2 ) )
    allocate( albg_org (       IALL, JALL, 2 ) )
    allocate( z0w_org  (       IALL, JALL    ) )
    allocate( skint_org(       IALL, JALL    ) )
    allocate( skinw_org(       IALL, JALL    ) )
    allocate( snowq_org(       IALL, JALL    ) )
    allocate( snowt_org(       IALL, JALL    ) )

    allocate( hfact(        IA, JA, itp_nh         ) )
    allocate( vfact( LKMAX, IA, JA, itp_nh, itp_nv ) )
    allocate( igrd (        IA, JA, itp_nh         ) )
    allocate( jgrd (        IA, JA, itp_nh         ) )
    allocate( kgrd ( LKMAX, IA, JA, itp_nh, itp_nv ) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputSCALE-Surface]'

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
       yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

       xs = PARENT_IMAX * (xloc-1) + 1
       xe = PARENT_IMAX * xloc
       ys = PARENT_JMAX * (yloc-1) + 1
       ye = PARENT_JMAX * yloc

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

       call FileRead( read3D(:,:,:), BASENAME_ORG, "LAND_WATER", 1, rank )
       do k = 1, KALL
         strg_org(k,xs:xe,ys:ye) = read3D(:,:,k)
       end do

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_TEMP",     1, rank )
       tw_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "LAND_SFC_TEMP",  1, rank )
       lst_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "TS_URB",         1, rank )
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

       call FileRead( read2D(:,:), BASENAME_ORG, "OCEAN_SFC_Z0",   1, rank )
       z0w_org(xs:xe,ys:ye) = read2D(:,:)

       call FileRead( read2D(:,:), BASENAME_ORG, "SFC_TEMP",       1, rank )
       skint_org(xs:xe,ys:ye) = read2D(:,:)

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

    ! tentative process
    skinw_org(:,:) = 0.0_RP
    snowq_org(:,:) = 0.0_RP
    snowt_org(:,:) = TEM00

    call latlonz_interporation_fact( hfact  (:,:,:),     & ! [OUT]
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

    ! interpolation
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = 1, LKMAX-1
       tg  (k,i,j) = tg_org  (kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                   + tg_org  (kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                   + tg_org  (kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                   + tg_org  (kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                   + tg_org  (kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                   + tg_org  (kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,2)
       strg(k,i,j) = strg_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                   + strg_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                   + strg_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                   + strg_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                   + strg_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                   + strg_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    enddo
    enddo
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       tg  (LKMAX,i,j) = tg  (LKMAX-1,i,j)
       strg(LKMAX,i,j) = strg(LKMAX-1,i,j)
    enddo
    enddo

    roff(:,:) = 0.0_RP ! not necessary
    qvef(:,:) = 0.0_RP ! not necessary

    do j = JS-1, JE+1
    do i = IS-1, IE+1
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

    do n = 1, 2 ! 1:LW, 2:SW
    do j = JS-1, JE+1
    do i = IS-1, IE+1
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

    deallocate( hfact )
    deallocate( vfact )
    deallocate( igrd  )
    deallocate( jgrd  )
    deallocate( kgrd  )

    return
  end subroutine InputSurfaceSCALE

  !-----------------------------------------------------------------------------
  subroutine InputSurfaceWRF( &
      tg,          & ! (out)
      strg,        & ! (out)
      roff,        & ! (out)
      qvef,        & ! (out)
      tw,          & ! (out)
      lst,         & ! (out)
      ust,         & ! (out)
      sst,         & ! (out)
      albw,        & ! (out)
      albg,        & ! (out)
      z0w,         & ! (out)
      skint,       & ! (out)
      skinw,       & ! (out)
      snowq,       & ! (out)
      snowt,       & ! (out)
      basename,    & ! (in)
      dims,        & ! (in)
      mdlid,       & ! (in)
      serial,      & ! (in)
      no_input     & ! (in)
      )
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
    character(LEN=*), intent( in) :: basename
    integer,          intent( in)  :: mdlid
    integer,          intent( in)  :: dims(:)
    logical,          intent( in)  :: serial
    logical,          intent( in)  :: no_input

    real(RP), allocatable :: org_3D(:,:,:)
    real(RP), allocatable :: org_4D(:,:,:,:)
    real(RP), allocatable :: lat_ORG(:,:,:)
    real(RP), allocatable :: lon_ORG(:,:,:)
    real(RP), allocatable :: zs_org(:,:,:,:)
    real(RP), allocatable :: dzs_org(:,:)

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
    integer :: n, k, i, j, iq, iqw
    integer :: tcount = 1

    logical :: do_read

    intrinsic shape
    !---------------------------------------------------------------------------
    d2r = pi / 180.0_RP

    allocate( lcz_3D(LKMAX,IA,JA) )
    allocate( hfact(IA,JA,itp_nh) )
    allocate( vfact(LKMAX,IA,JA,itp_nh,itp_nv) )
    allocate( igrd(IA,JA,itp_nh) )
    allocate( jgrd(IA,JA,itp_nh) )
    allocate( kgrd(LKMAX,IA,JA,itp_nh,itp_nv) )

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

    if( IO_L ) write(IO_FID_LOG,*)
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
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,lcz_3D,lat,lon, &
                                      zs_org,lat_org,lon_org,dims(7),dims(2),dims(3),1,landgrid=.true. )

    ! soil temperature [K]
    if( do_read ) then
       call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                      BASENAME, "TSLB",  1, tcount, myrank, mdlid, single=.true., landgrid=.true. )
       org_4D(:,:,:,:) = dummy_4D(:,:,:,:)
    endif
    if( serial ) call COMM_bcast( org_4D(:,:,:,:), dims(7),dims(2),dims(3),tcount )
    n = 1
    do j = JS-1, JE+1
    do i = IS-1, IE+1
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
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       roff(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                  + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                  + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ! soil liquid water m3 m-3] (no wrfout-default)
    ! --- not available; appropriate convert method is not found
    strg(:,:,:) = 0.0_DP

    ! accumulated surface evaporation [kg m-2]
    ! --- not available; convert method is not found
    qvef(:,:) = 0.0_DP

    ! sea surface temperature [K]
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "SST",  1, tcount, myrank, mdlid, single=.true. )
       org_3D(:,:,:) = dummy_3D(:,:,:)
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    n = 1
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       tw(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ! surface skin temperature [K]
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "TSK",  1, tcount, myrank, mdlid, single=.true. )
       org_3D(:,:,:) = dummy_3D(:,:,:)
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    n = 1
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       lst(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                 + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                 + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ust(:,:) = lst(:,:)
    sst(:,:) = lst(:,:)

    ! ALBEDO [-]
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "ALBEDO",  1, tcount, myrank, mdlid, single=.true. )
       org_3D(:,:,:) = dummy_3D(:,:,:)
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    n = 1
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       albw(i,j,I_SW) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                       + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                       + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    albg(:,:,I_SW) = albw(:,:,I_SW)

    ! SURFACE EMISSIVITY [-] (no wrfout-default)
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
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       albw(i,j,I_LW) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                       + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                       + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    albg(:,:,I_LW) = albw(:,:,I_LW)

    ! TIME-VARYING ROUGHNESS LENGTH [m] (no wrfout-default)
    if ( no_input ) then
       z0w(:,:) = 0.1_DP
    else
       if( do_read ) then
          call ExternalFileRead( dummy_3D(:,:,:),                             &
                         BASENAME, "ZNT",  1, tcount, myrank, mdlid, single=.true. )
          org_3D(:,:,:) = dummy_3D(:,:,:)
       endif
       if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
       n = 1
       do j = JS-1, JE+1
       do i = IS-1, IE+1
          z0w(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                    + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                    + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
       enddo
       enddo
    endif

    !tentative approach for skin
    skint(:,:) = lst(:,:)
    skinw(:,:) = 0.0_DP

    ! SNOW WATER EQUIVALENT [kg m-2] (no wrfout-default)
    if( do_read ) then
       call ExternalFileRead( dummy_3D(:,:,:),                             &
                      BASENAME, "SNOW",  1, tcount, myrank, mdlid, single=.true. )
       org_3D(:,:,:) = dummy_3D(:,:,:)
    endif
    if( serial ) call COMM_bcast( org_3D(:,:,:), dims(2),dims(3),tcount )
    n = 1
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       snowq(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                   + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                   + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ! AVERAGE SNOW TEMPERATURE [C] (no wrfout-default)
    if ( no_input ) then
       snowt(:,:) = lst(:,:)
    else
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
       do j = JS-1, JE+1
       do i = IS-1, IE+1
          snowt(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                      + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                      + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
       enddo
       enddo
    endif

    deallocate( dummy_2D )
    deallocate( dummy_3D )
    deallocate( dummy_4D )
    deallocate( org_3D )
    deallocate( org_4D )

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
  subroutine LatLonZ_Interpolation_Linear( &
      dens,          & ! (out)
      w,             & ! (out)
      u,             & ! (out)
      v,             & ! (out)
      pott,          & ! (out)
      qtrc,          & ! (out)
      psfc,          & ! (out)
      th2m,          & ! (out)
      dens_org,      & ! (in)
      w_org,         & ! (in)
      u_org,         & ! (in)
      v_org,         & ! (in)
      pott_org,      & ! (in)
      qtrc_org,      & ! (in)
      psfc_org,      & ! (in)
      th2m_org,      & ! (in)
      lat_org,       & ! (in)
      lon_org,       & ! (in)
      latu_org,      & ! (in)
      lonu_org,      & ! (in)
      latv_org,      & ! (in)
      lonv_org,      & ! (in)
      geof_org,      & ! (in)
      geoh_org,      & ! (in)
      dims,          & ! (in)
      ts,            & ! (in)
      te             & ! (in)
      )
    implicit none

    real(RP), intent(out)  :: dens(:,:,:,:)
    real(RP), intent(out)  :: w(:,:,:,:)
    real(RP), intent(out)  :: u(:,:,:,:)
    real(RP), intent(out)  :: v(:,:,:,:)
    real(RP), intent(out)  :: pott(:,:,:,:)
    real(RP), intent(out)  :: qtrc(:,:,:,:,:)
    real(RP), intent(out)  :: psfc(:,:,:)
    real(RP), intent(out)  :: th2m(:,:,:)
    real(RP), intent( in)  :: dens_org(:,:,:,:)
    real(RP), intent( in)  :: w_org(:,:,:,:)
    real(RP), intent( in)  :: u_org(:,:,:,:)
    real(RP), intent( in)  :: v_org(:,:,:,:)
    real(RP), intent( in)  :: pott_org(:,:,:,:)
    real(RP), intent( in)  :: qtrc_org(:,:,:,:,:)
    real(RP), intent( in)  :: psfc_org(:,:,:)
    real(RP), intent( in)  :: th2m_org(:,:,:)
    real(RP), intent( in)  :: lat_org(:,:,:)
    real(RP), intent( in)  :: lon_org(:,:,:)
    real(RP), intent( in)  :: latu_org(:,:,:)
    real(RP), intent( in)  :: lonu_org(:,:,:)
    real(RP), intent( in)  :: latv_org(:,:,:)
    real(RP), intent( in)  :: lonv_org(:,:,:)
    real(RP), intent( in)  :: geof_org(:,:,:,:)
    real(RP), intent( in)  :: geoh_org(:,:,:,:)
    integer,  intent( in)  :: dims(:)
    integer,  intent( in)  :: ts
    integer,  intent( in)  :: te

    real(RP) :: hfact(IA,JA,itp_nh)
    real(RP) :: vfact(KA,IA,JA,itp_nh,itp_nv)
    integer :: igrd(IA,JA,itp_nh)
    integer :: jgrd(IA,JA,itp_nh)
    integer :: kgrd(KA,IA,JA,itp_nh,itp_nv)

    integer :: n, k, i, j, iq

    intrinsic shape
    !---------------------------------------------------------------------------

    do n=ts, te !--- time loop

    ! for scalar points
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,cz,lat,lon, &
                                      geoh_org,lat_org,lon_org,dims(1),dims(2),dims(3),n )

    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS-1, KE+1
       dens(k,i,j,n) =  dens_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                      + dens_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                      + dens_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                      + dens_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                      + dens_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                      + dens_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)

       pott(k,i,j,n) =  pott_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                      + pott_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                      + pott_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                      + pott_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                      + pott_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                      + pott_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)

       do iq = 1, QA
          qtrc(k,i,j,n,iq) =  qtrc_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n,iq) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                            + qtrc_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n,iq) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                            + qtrc_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n,iq) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                            + qtrc_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n,iq) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                            + qtrc_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n,iq) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                            + qtrc_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n,iq) * hfact(i,j,3) * vfact(k,i,j,3,2)
       enddo
    enddo
    enddo
    enddo

    do j = JS-1, JE+1
    do i = IS-1, IE+1
       psfc(i,j,n) =    psfc_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                      + psfc_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                      + psfc_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)

       th2m(i,j,n)  =   th2m_org(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                      + th2m_org(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                      + th2m_org(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    ! for vector (u) points
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,cz,lat,lon_u, &
                                      geoh_org,latu_org,lonu_org,dims(1),dims(2),dims(3),n )
    !                                 dims(5) wasn't used to keep consistency with geoh_org
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS-1, KE+1
       u(k,i,j,n) =  u_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                   + u_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                   + u_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                   + u_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                   + u_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                   + u_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    enddo
    enddo

    ! for vector (v) points
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,cz,lat_v,lon, &
                                      geoh_org,latv_org,lonv_org,dims(1),dims(2),dims(3),n )
    !                                 dims(6) wasn't used to keep consistency with geoh_org
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS-1, KE+1
       v(k,i,j,n) =  v_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                   + v_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                   + v_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                   + v_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                   + v_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                   + v_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    enddo
    enddo

    ! for vector (w) points
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,fz,lat,lon, &
                                      geof_org,lat_org,lon_org,dims(4),dims(2),dims(3),n )
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS-1, KE+1
       w(k,i,j,n) =  w_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                   + w_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                   + w_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                   + w_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                   + w_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                   + w_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    enddo
    enddo

    enddo !--- time loop

    return
  end subroutine LatLonZ_Interpolation_Linear

  !-----------------------------------------------------------------------------
  subroutine latlonz_interporation_fact( &
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

    do j = JS-1, JE+1
    do i = IS-1, IE+1
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
  end subroutine latlonz_interporation_fact

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

    qvars(:,:,:,:,I_NC) = qvars(:,:,:,:,I_QC) / ( (piov6*RHOw) * Dc**b )
    qvars(:,:,:,:,I_NR) = qvars(:,:,:,:,I_QR) / ( (piov6*RHOw) * Dr**b )
    qvars(:,:,:,:,I_NI) = qvars(:,:,:,:,I_QI) / ( (piov6*RHOf) * Di**b )
    qvars(:,:,:,:,I_NS) = qvars(:,:,:,:,I_QS) / ( (piov6*RHOf) * Ds**b )
    qvars(:,:,:,:,I_NG) = qvars(:,:,:,:,I_QG) / ( (piov6*RHOg) * Dg**b )

    return
  end subroutine diagnose_number_concentration

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
