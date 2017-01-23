!-------------------------------------------------------------------------------
!> module REAL input
!!
!! @par Description
!!          read data from file for real atmospheric simulations
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
  use scale_land_grid_index
  use scale_urban_grid_index
  use scale_grid_real, only: &
     LON  => REAL_LON,  &
     LAT  => REAL_LAT,  &
     LONX => REAL_LONX, &
     LATY => REAL_LATY, &
     CZ   => REAL_CZ,   &
     FZ   => REAL_FZ
  use scale_grid_nest, only: &
     NEST_INTERP_LEVEL
  use scale_index
  use scale_tracer
  use scale_process, only: &
     PRC_IsMaster, &
     PRC_MPIstop
  use scale_external_io, only: &
     iSCALE, &
     iWRFARW, &
     iNICAM, &
     iGrADS
  use scale_comm, only: &
     COMM_bcast
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: REALINPUT_Atmos
  public :: REALINPUT_Surface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ParentAtomInput
  private :: ParentAtomBoundary
  private :: ParentSurfaceSetup
  private :: ParentSurfaceInput
  private :: ParentOceanBoundary
  private :: interp_OceanLand_data

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: I_intrp_off  = 0
  integer, parameter :: I_intrp_mask = 1
  integer, parameter :: I_intrp_fill = 2

  integer, parameter :: cosin = 1
  integer, parameter :: sine  = 2

  real(RP), allocatable :: lon_org (:,:)
  real(RP), allocatable :: lat_org (:,:)
  real(RP), allocatable :: cz_org(:,:,:)

  real(RP), allocatable :: dens_org(:,:,:)
  real(RP), allocatable :: qtrc_org(:,:,:,:)

  real(RP), allocatable :: velz_org(:,:,:)
  real(RP), allocatable :: velx_org(:,:,:)
  real(RP), allocatable :: vely_org(:,:,:)
  real(RP), allocatable :: pott_org(:,:,:)
  real(RP), allocatable :: temp_org(:,:,:)
  real(RP), allocatable :: pres_org(:,:,:)

  real(RP), allocatable :: hfact(:,:,:)
  real(RP), allocatable :: vfact(:,:,:,:,:)
  integer,  allocatable :: igrd (:,:,:)
  integer,  allocatable :: jgrd (:,:,:)
  integer,  allocatable :: kgrd (:,:,:,:,:)
  integer,  allocatable :: ncopy(:,:,:)

  real(RP), allocatable :: tw_org(:,:)
  real(RP), allocatable :: sst_org(:,:)
  real(RP), allocatable :: albw_org(:,:,:)
  real(RP), allocatable :: olon_org(:,:)
  real(RP), allocatable :: olat_org(:,:)
  real(RP), allocatable :: omask_org(:,:)

  integer, private   :: itp_nh = 4
  integer, private   :: itp_nv = 2

  integer, private   :: io_fid_grads_nml  = -1
  integer, private   :: io_fid_grads_data = -1

  logical, private   :: do_read_atom
  logical, private   :: do_read_land
  logical, private   :: do_read_ocean
  logical, private   :: rotate
  logical, private   :: use_waterratio
  logical, private   :: update_coord
  logical, private   :: use_temp
  logical, private   :: serial
  logical, private   :: serial_land
  logical, private   :: serial_ocean
  logical, private   :: first = .true.

  integer, private   :: i_intrp_land_temp
  integer, private   :: i_intrp_land_water
  integer, private   :: i_intrp_land_sfc_temp
  integer, private   :: i_intrp_ocean_temp
  integer, private   :: i_intrp_ocean_sfc_temp


  ! replace missing value
  real(RP), private, parameter :: maskval_tg   = 298.0_RP !> mask value 298K
  real(RP), private, parameter :: maskval_strg = 0.02_RP  !> mask value 0.02
                                ! default value 0.02: set as value of forest at 40% of evapolation rate.
                                ! forest is considered as a typical landuse over Japan area.


  ! for namelist
  integer                  :: NUMBER_OF_FILES     = 1
  integer                  :: NUMBER_OF_TSTEPS    = 1    ! num of time steps in one file
  integer                  :: NUMBER_OF_SKIP_TSTEPS  = 0 ! num of skipped first several data

  character(len=H_LONG)    :: FILETYPE_ORG        = ''
  character(len=H_LONG)    :: BASENAME_ORG        = ''
  logical                  :: BASENAME_ADD_NUM    = .false.
  character(len=H_LONG)    :: BASENAME_BOUNDARY   = 'boundary_atmos'
  character(len=H_LONG)    :: BOUNDARY_TITLE      = 'SCALE-RM BOUNDARY CONDITION for REAL CASE'
  real(RP)                 :: BOUNDARY_UPDATE_DT  = 0.0_RP    ! inteval time of boudary data update [s]

  integer                  :: PARENT_MP_TYPE      = 6         ! microphysics type of the parent model (number of classes)
                                                                ! 0: dry, 3:3-class, 5:5-class, 6:6-class, >6:double moment
  logical                  :: SERIAL_PROC_READ    = .true.    ! read by one MPI process and broadcast
  ! only for SCALE boundary
  logical                  :: USE_FILE_DENSITY    = .false.   ! use density data from files
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine REALINPUT_atmos( &
       flg_intrp )
    use mod_atmos_vars, only: &
         DENS, &
         MOMZ, &
         MOMX, &
         MOMY, &
         RHOT, &
         QTRC
    use mod_atmos_admin, only: &
         ATMOS_PHY_MP_TYPE
    implicit none

    logical, intent(in)  :: flg_intrp ! flag for interpolation of SBM(S10) from outer bulk-MP model


    NAMELIST / PARAM_MKINIT_REAL_ATMOS / &
         NUMBER_OF_FILES,        &
         NUMBER_OF_TSTEPS,       &
         NUMBER_OF_SKIP_TSTEPS,  &
         FILETYPE_ORG,           &
         BASENAME_ORG,           &
         BASENAME_ADD_NUM,       &
         BASENAME_BOUNDARY,      &
         BOUNDARY_TITLE,         &
         BOUNDARY_UPDATE_DT,     &
         PARENT_MP_TYPE,         &
         SERIAL_PROC_READ,       &
         USE_FILE_DENSITY

    character(len=H_LONG) :: BASENAME_ATMOS    = ''
    character(len=H_LONG) :: BASENAME_WITHNUM  = ''
    character(len=5)      :: NUM               = ''

    ! atmos
    real(RP), allocatable :: DENS_ORG(:,:,:,:)
    real(RP), allocatable :: MOMZ_ORG(:,:,:,:)
    real(RP), allocatable :: MOMX_ORG(:,:,:,:)
    real(RP), allocatable :: MOMY_ORG(:,:,:,:)
    real(RP), allocatable :: RHOT_ORG(:,:,:,:)
    real(RP), allocatable :: QTRC_ORG(:,:,:,:,:)

    integer :: mdlid
    integer :: dims(6) ! dims 1-3: normal, 4-6: staggerd

    integer :: totaltimesteps = 1
    integer :: timelen
    integer :: skip_steps
    integer :: ierr
    logical :: flg_bin          ! flag for SBM(S10) is used or not 0-> not used, 1-> used

    integer :: k, i, j, iq, n, ns, ne
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[RealCaseAtmos]/Categ[MKINIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_ATMOS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_ATMOS. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_MKINIT_REAL_ATMOS)

    if( ATMOS_PHY_MP_TYPE == 'SUZUKI10' ) then
       flg_bin = .true.
    else
       flg_bin = .false.
    endif

    if ( FILETYPE_ORG == "GrADS" ) then
       BASENAME_WITHNUM = BASENAME_ORG ! namelist file name
       BASENAME_ATMOS = ""
    else
       if ( NUMBER_OF_FILES > 1 .or. BASENAME_ADD_NUM ) then
          BASENAME_WITHNUM = trim(BASENAME_ORG)//"_00000"
       else
          BASENAME_WITHNUM = trim(BASENAME_ORG)
       end if
       BASENAME_ATMOS = BASENAME_ORG
    end if

    call ParentAtomSetup( dims(:), timelen, mdlid, & ![OUT]
                          BASENAME_WITHNUM,        & ![IN]
                          FILETYPE_ORG,            & ![IN]
                          USE_FILE_DENSITY,        & ![IN]
                          SERIAL_PROC_READ         ) ![IN]

    if ( BOUNDARY_UPDATE_DT <= 0.0_RP ) then
       write(*,*) 'xxx BOUNDARY_UPDATE_DT is necessary in real case preprocess'
       call PRC_MPIstop
    endif

    if ( timelen > 0 ) then
       NUMBER_OF_TSTEPS = timelen ! read from file
    endif

    totaltimesteps = NUMBER_OF_FILES * NUMBER_OF_TSTEPS

    allocate( dens_org(KA,IA,JA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )
    allocate( momz_org(KA,IA,JA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )
    allocate( momx_org(KA,IA,JA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )
    allocate( momy_org(KA,IA,JA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )
    allocate( rhot_org(KA,IA,JA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )
    allocate( qtrc_org(KA,IA,JA,QA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )

    !--- read external file
    do n = 1, NUMBER_OF_FILES

       if ( NUMBER_OF_FILES > 1 .or. BASENAME_ADD_NUM ) then
          write(NUM,'(I5.5)') n-1
          BASENAME_WITHNUM = trim(BASENAME_ATMOS)//"_"//NUM
       else
          BASENAME_WITHNUM = trim(BASENAME_ATMOS)
       end if

       if( IO_L ) write(IO_FID_LOG,*) ' '
       if( IO_L ) write(IO_FID_LOG,*) '+++ Target File Name: ',trim(BASENAME_WITHNUM)
       if( IO_L ) write(IO_FID_LOG,*) '    Time Steps in One File: ', NUMBER_OF_TSTEPS

       ns = NUMBER_OF_TSTEPS * (n - 1) + 1
       ne = ns + (NUMBER_OF_TSTEPS - 1)

       if ( ne <= NUMBER_OF_SKIP_TSTEPS ) then
          if( IO_L ) write(IO_FID_LOG,*) '    SKIP'
          cycle
       end if

       skip_steps = max(NUMBER_OF_SKIP_TSTEPS - ns + 1, 0)
       ns = max(ns, NUMBER_OF_SKIP_TSTEPS+1)

       ! read all prepared data
       call ParentAtomInput( DENS_org(:,:,:,ns:ne),   &
                             MOMZ_org(:,:,:,ns:ne),   &
                             MOMX_org(:,:,:,ns:ne),   &
                             MOMY_org(:,:,:,ns:ne),   &
                             RHOT_org(:,:,:,ns:ne),   &
                             QTRC_org(:,:,:,:,ns:ne), &
                             BASENAME_WITHNUM,        &
                             dims(:),                 &
                             mdlid,                   &
                             flg_bin, flg_intrp,      &
                             PARENT_MP_TYPE,          &
                             NUMBER_OF_TSTEPS,        &
                             skip_steps               )
    enddo

    !--- input initial data
    ns = NUMBER_OF_SKIP_TSTEPS + 1  ! skip first several data
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       DENS(k,i,j) = DENS_ORG(k,i,j,ns)
       MOMZ(k,i,j) = MOMZ_ORG(k,i,j,ns)
       MOMX(k,i,j) = MOMX_ORG(k,i,j,ns)
       MOMY(k,i,j) = MOMY_ORG(k,i,j,ns)
       RHOT(k,i,j) = RHOT_ORG(k,i,j,ns)

       do iq = 1, QA
          QTRC(k,i,j,iq) = QTRC_ORG(k,i,j,iq,ns)
       enddo
    enddo
    enddo
    enddo

    !--- output boundary data
    totaltimesteps = totaltimesteps - NUMBER_OF_SKIP_TSTEPS ! skip first several data
    call ParentAtomBoundary( DENS_ORG(:,:,:,ns:ne),   &
                             MOMZ_ORG(:,:,:,ns:ne),   &
                             MOMX_ORG(:,:,:,ns:ne),   &
                             MOMY_ORG(:,:,:,ns:ne),   &
                             RHOT_ORG(:,:,:,ns:ne),   &
                             QTRC_ORG(:,:,:,:,ns:ne), &
                             totaltimesteps,          &
                             BOUNDARY_UPDATE_DT,      &
                             BASENAME_BOUNDARY,       &
                             BOUNDARY_TITLE           )

    deallocate( dens_org )
    deallocate( momz_org )
    deallocate( momx_org )
    deallocate( momy_org )
    deallocate( rhot_org )
    deallocate( qtrc_org )

    return
  end subroutine REALINPUT_atmos

  !-----------------------------------------------------------------------------
  subroutine REALINPUT_surface
    use scale_const, only: &
         I_SW => CONST_I_SW, &
         I_LW => CONST_I_LW
    use mod_atmos_vars, only: &
         DENS, &
         MOMZ, &
         MOMX, &
         MOMY, &
         RHOT, &
         QTRC
    use mod_land_vars, only: &
         LAND_TEMP, &
         LAND_WATER, &
         LAND_SFC_TEMP, &
         LAND_SFC_albedo
    use mod_urban_vars, only: &
         URBAN_SFC_TEMP, &
         URBAN_SFC_albedo, &
         URBAN_TC, &
         URBAN_QC, &
         URBAN_UC, &
         URBAN_TR, &
         URBAN_TB, &
         URBAN_TG, &
         URBAN_TRL, &
         URBAN_TBL, &
         URBAN_TGL, &
         URBAN_RAINR, &
         URBAN_RAINB, &
         URBAN_RAING, &
         URBAN_ROFF
    use mod_ocean_vars, only: &
         OCEAN_TEMP, &
         OCEAN_SFC_TEMP, &
         OCEAN_SFC_albedo, &
         OCEAN_SFC_Z0M, &
         OCEAN_SFC_Z0H, &
         OCEAN_SFC_Z0E
    use mod_atmos_phy_sf_vars, only: &
         ATMOS_PHY_SF_SFC_TEMP, &
         ATMOS_PHY_SF_SFC_albedo, &
         ATMOS_PHY_SF_SFC_Z0M, &
         ATMOS_PHY_SF_SFC_Z0H, &
         ATMOS_PHY_SF_SFC_Z0E
    use scale_landuse, only: &
         fact_ocean  => LANDUSE_fact_ocean, &
         fact_land   => LANDUSE_fact_land, &
         fact_urban  => LANDUSE_fact_urban
    implicit none

    logical                  :: USE_FILE_LANDWATER   = .true.  ! use land water data from files
    real(RP)                 :: INIT_LANDWATER_RATIO = 0.5_RP  ! Ratio of land water to storage is constant, if USE_FILE_LANDWATER is ".false."
    character(len=H_SHORT)   :: INTRP_LAND_TEMP      = 'off'
    character(len=H_SHORT)   :: INTRP_LAND_WATER     = 'off'
    character(len=H_SHORT)   :: INTRP_LAND_SFC_TEMP  = 'off'
    character(len=H_SHORT)   :: INTRP_OCEAN_TEMP     = 'off'
    character(len=H_SHORT)   :: INTRP_OCEAN_SFC_TEMP = 'off'
    integer                  :: INTRP_ITER_MAX       = 20
    character(len=H_SHORT)   :: SOILWATER_DS2VC      = 'limit'
    logical                  :: soilwater_DS2VC_flag           ! true: 'critical', false: 'limit'
    logical                  :: elevation_collection = .true.

    NAMELIST / PARAM_MKINIT_REAL_LAND / &
         NUMBER_OF_FILES,        &
         NUMBER_OF_TSTEPS,       &
         NUMBER_OF_SKIP_TSTEPS,  &
         FILETYPE_ORG,           &
         BASENAME_ORG,           &
         BASENAME_ADD_NUM,       &
         USE_FILE_LANDWATER,     &
         INIT_LANDWATER_RATIO,   &
         INTRP_LAND_TEMP,        &
         INTRP_LAND_WATER,       &
         INTRP_LAND_SFC_TEMP,    &
         INTRP_ITER_MAX,         &
         SOILWATER_DS2VC,        &
         ELEVATION_COLLECTION,   &
         SERIAL_PROC_READ

    NAMELIST / PARAM_MKINIT_REAL_OCEAN / &
         NUMBER_OF_FILES,        &
         NUMBER_OF_TSTEPS,       &
         NUMBER_OF_SKIP_TSTEPS,  &
         FILETYPE_ORG,           &
         BASENAME_ORG,           &
         BASENAME_ADD_NUM,       &
         BASENAME_BOUNDARY,      &
         BOUNDARY_TITLE,         &
         BOUNDARY_UPDATE_DT,     &
         INTRP_OCEAN_TEMP,       &
         INTRP_OCEAN_SFC_TEMP,   &
         INTRP_ITER_MAX,         &
         SERIAL_PROC_READ

    character(len=H_LONG) :: FILETYPE_LAND
    character(len=H_LONG) :: FILETYPE_OCEAN
    character(len=H_LONG) :: BASENAME_LAND
    character(len=H_LONG) :: BASENAME_OCEAN
    character(len=H_LONG) :: BASENAME_WITHNUM  = ''
    character(len=5)      :: NUM               = ''
    logical               :: SERIAL_PROC_READ_land
    logical               :: SERIAL_PROC_READ_ocean

    ! land
    real(RP) :: LAND_TEMP_ORG(LKMAX,IA,JA)
    real(RP) :: LAND_WATER_ORG(LKMAX,IA,JA)
    real(RP) :: LAND_SFC_TEMP_ORG(IA,JA)
    real(RP) :: LAND_SFC_albedo_ORG(IA,JA,2)

    ! urban
    real(RP) :: URBAN_TC_ORG(IA,JA)
    real(RP) :: URBAN_QC_ORG(IA,JA)
    real(RP) :: URBAN_UC_ORG(IA,JA)
    real(RP) :: URBAN_SFC_TEMP_ORG(IA,JA)
    real(RP) :: URBAN_SFC_albedo_ORG(IA,JA,2)

    ! ocean
    real(RP), allocatable :: OCEAN_TEMP_ORG(:,:,:)
    real(RP), allocatable :: OCEAN_SFC_TEMP_ORG(:,:,:)
    real(RP), allocatable :: OCEAN_SFC_albedo_ORG(:,:,:,:)
    real(RP), allocatable :: OCEAN_SFC_Z0_ORG(:,:,:)

    integer :: mdlid_land, mdlid_ocean
    integer :: ldims(3), odims(2)

    integer :: totaltimesteps = 1
    integer :: timelen
    integer :: skip_steps
    integer :: lit
    integer :: lfn
    integer :: ierr

    integer :: k, i, j, n, ns, ne
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[RealCaseSurface]/Categ[MKINIT]'


    ! LAND/URBAN

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_LAND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_LAND. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_MKINIT_REAL_LAND)

    FILETYPE_LAND = FILETYPE_ORG

    lfn = NUMBER_OF_SKIP_TSTEPS / NUMBER_OF_TSTEPS
    if ( FILETYPE_LAND .ne. "GrADS" .and. ( NUMBER_OF_FILES > 1 .or. BASENAME_ADD_NUM ) ) then
       write(NUM,'(I5.5)') lfn
       BASENAME_LAND = trim(BASENAME_ORG)//"_"//NUM
    else
       BASENAME_LAND = trim(BASENAME_ORG)
    end if

    serial_land = SERIAL_PROC_READ

    lit = mod(NUMBER_OF_SKIP_TSTEPS,NUMBER_OF_TSTEPS)+1

    !--- read external file
    if( IO_L ) write(IO_FID_LOG,*) ' '
    if( IO_L ) write(IO_FID_LOG,*) '+++ Target File Name (Land): ',trim(BASENAME_LAND)
    if( IO_L ) write(IO_FID_LOG,*) '    Time Steps in One File: ', NUMBER_OF_TSTEPS
    if( IO_L ) write(IO_FID_LOG,*) '    Time Step to read: ', lit



    ! OCEAN

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_OCEAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_OCEAN. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_MKINIT_REAL_OCEAN)

    FILETYPE_OCEAN = FILETYPE_ORG

    if ( FILETYPE_OCEAN .ne. "GrADS" .and. ( NUMBER_OF_FILES > 1 .or. BASENAME_ADD_NUM ) ) then
       BASENAME_OCEAN = trim(BASENAME_ORG)//"_00000"
    else
       BASENAME_OCEAN = trim(BASENAME_ORG)
    end if

    select case( SOILWATER_DS2VC )
    case( 'critical' )
       SOILWATER_DS2VC_flag = .true.
    case('limit' )
       SOILWATER_DS2VC_flag = .false.
    case default
      write(*,*) 'xxx Unsupported SOILWATER_DS2CV TYPE:', trim(SOILWATER_DS2VC)
      call PRC_MPIstop
    end select

    serial_ocean = SERIAL_PROC_READ


    call ParentSurfaceSetup( ldims, odims,           & ![OUT]
                             mdlid_land,             & ![OUT]
                             mdlid_ocean,            & ![OUT]
                             timelen,                & ![OUT]
                             BASENAME_LAND,          & ![IN]
                             BASENAME_OCEAN,         & ![IN]
                             FILETYPE_LAND,          & ![IN]
                             FILETYPE_OCEAN,         & ![IN]
                             USE_FILE_LANDWATER,     & ![IN]
                             intrp_land_temp,        & ![IN]
                             intrp_land_water,       & ![IN]
                             intrp_land_sfc_temp,    & ![IN]
                             intrp_ocean_temp,       & ![IN]
                             intrp_ocean_sfc_temp    ) ![IN]

    if ( timelen > 0 ) then
       NUMBER_OF_TSTEPS = timelen ! read from file
    endif

    totaltimesteps = NUMBER_OF_FILES * NUMBER_OF_TSTEPS

    allocate( OCEAN_TEMP_ORG      (IA,JA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )
    allocate( OCEAN_SFC_TEMP_ORG  (IA,JA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )
    allocate( OCEAN_SFC_albedo_ORG(IA,JA,2,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )
    allocate( OCEAN_SFC_Z0_ORG    (IA,JA,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps   ) )

    if ( mdlid_land == iGrADS .and. ( NUMBER_OF_FILES > 1 .or. BASENAME_ADD_NUM ) ) then
       write(NUM,'(I5.5)') lfn
       BASENAME_LAND = "_"//NUM
    end if

    if ( mdlid_ocean == iGrADS ) then
       BASENAME_ORG = ""
    end if

    !--- read external file
    do n = 1, NUMBER_OF_FILES

       if ( NUMBER_OF_FILES > 1 .or. BASENAME_ADD_NUM ) then
          write(NUM,'(I5.5)') n-1
          BASENAME_OCEAN = trim(BASENAME_ORG)//"_"//NUM
       else
          BASENAME_OCEAN = trim(BASENAME_ORG)
       end if

       if( IO_L ) write(IO_FID_LOG,*) ' '
       if( IO_L ) write(IO_FID_LOG,*) '+++ Target File Name (Ocean): ', trim(BASENAME_OCEAN)
       if( IO_L ) write(IO_FID_LOG,*) '    Time Steps in One File: ', NUMBER_OF_TSTEPS

       ns = NUMBER_OF_TSTEPS * (n - 1) + 1
       ne = ns + (NUMBER_OF_TSTEPS - 1)

       if ( ne <= NUMBER_OF_SKIP_TSTEPS ) then
          if( IO_L ) write(IO_FID_LOG,*) '    SKIP'
          cycle
       end if

       skip_steps = max(NUMBER_OF_SKIP_TSTEPS - ns + 1, 0)
       ns = max(ns, NUMBER_OF_SKIP_TSTEPS+1)

       ! read all prepared data
       call ParentSurfaceInput( LAND_TEMP_org,        &
                                LAND_WATER_org,       &
                                LAND_SFC_TEMP_org,    &
                                LAND_SFC_albedo_org,  &
                                URBAN_TC_org,         &
                                URBAN_QC_org,         &
                                URBAN_UC_org,         &
                                URBAN_SFC_TEMP_org,   &
                                URBAN_SFC_albedo_org, &
                                OCEAN_TEMP_org(:,:,ns:ne),         &
                                OCEAN_SFC_TEMP_org(:,:,ns:ne),     &
                                OCEAN_SFC_albedo_org(:,:,:,ns:ne), &
                                OCEAN_SFC_Z0_org(:,:,ns:ne),       &
                                DENS, &
                                MOMZ, &
                                MOMX, &
                                MOMY, &
                                RHOT, &
                                QTRC, &
                                BASENAME_LAND,           &
                                BASENAME_OCEAN,          &
                                mdlid_land, mdlid_ocean, &
                                ldims, odims,            &
                                USE_FILE_LANDWATER,      &
                                INIT_LANDWATER_RATIO,    &
                                INTRP_ITER_MAX,          &
                                SOILWATER_DS2VC_flag,    &
                                elevation_collection,    &
                                NUMBER_OF_TSTEPS,        &
                                skip_steps, lit          )

    enddo


    !--- input initial data
    do j = 1, JA
    do i = 1, IA
       LAND_SFC_TEMP(i,j) = LAND_SFC_TEMP_org(i,j)
       LAND_SFC_albedo(i,j,I_LW) = LAND_SFC_albedo_org(i,j,I_LW)
       LAND_SFC_albedo(i,j,I_SW) = LAND_SFC_albedo_org(i,j,I_SW)
       do k = 1, LKMAX
          LAND_TEMP(k,i,j) = LAND_TEMP_org(k,i,j)
          LAND_WATER(k,i,j) = LAND_WATER_org(k,i,j)
       end do

       URBAN_SFC_TEMP(i,j) = URBAN_SFC_TEMP_org(i,j)
       URBAN_SFC_albedo(i,j,I_LW) = URBAN_SFC_albedo_org(i,j,I_LW)
       URBAN_SFC_albedo(i,j,I_SW) = URBAN_SFC_albedo_org(i,j,I_SW)
       do k = UKS, UKE
          URBAN_TRL(k,i,j) = URBAN_SFC_TEMP_org(i,j)
          URBAN_TBL(k,i,j) = URBAN_SFC_TEMP_org(i,j)
          URBAN_TGL(k,i,j) = URBAN_SFC_TEMP_org(i,j)
       end do
       URBAN_TC(i,j) = URBAN_TC_org(i,j)
       URBAN_QC(i,j) = URBAN_QC_org(i,j)
       URBAN_UC(i,j) = URBAN_UC_org(i,j)
       URBAN_TR(i,j) = URBAN_SFC_TEMP_org(i,j)
       URBAN_TB(i,j) = URBAN_SFC_TEMP_org(i,j)
       URBAN_TG(i,j) = URBAN_SFC_TEMP_org(i,j)
       URBAN_RAINR(i,j) = 0.0_RP
       URBAN_RAINB(i,j) = 0.0_RP
       URBAN_RAING(i,j) = 0.0_RP
       URBAN_ROFF (i,j) = 0.0_RP
    enddo
    enddo


    ns = NUMBER_OF_SKIP_TSTEPS + 1  ! skip first several data
    do j = 1, JA
    do i = 1, IA
       OCEAN_TEMP(i,j)     = OCEAN_TEMP_ORG(i,j,ns)
       OCEAN_SFC_TEMP(i,j) = OCEAN_SFC_TEMP_ORG(i,j,ns)
       OCEAN_SFC_albedo(i,j,I_LW) = OCEAN_SFC_albedo_ORG(i,j,I_LW,ns)
       OCEAN_SFC_albedo(i,j,I_SW) = OCEAN_SFC_albedo_ORG(i,j,I_SW,ns)
       OCEAN_SFC_Z0M(i,j) = OCEAN_SFC_Z0_ORG(i,j,ns)
       OCEAN_SFC_Z0H(i,j) = OCEAN_SFC_Z0_ORG(i,j,ns)
       OCEAN_SFC_Z0E(i,j) = OCEAN_SFC_Z0_ORG(i,j,ns)
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       ATMOS_PHY_SF_SFC_TEMP  (i,j)      = fact_ocean(i,j) * OCEAN_SFC_TEMP(i,j) &
                                         + fact_land (i,j) * LAND_SFC_TEMP (i,j) &
                                         + fact_urban(i,j) * URBAN_SFC_TEMP(i,j)
       ATMOS_PHY_SF_SFC_albedo(i,j,I_LW) = fact_ocean(i,j) * OCEAN_SFC_albedo(i,j,I_LW) &
                                         + fact_land (i,j) * LAND_SFC_albedo(i,j,I_LW) &
                                         + fact_urban(i,j) * URBAN_SFC_albedo(i,j,I_LW)
       ATMOS_PHY_SF_SFC_albedo(i,j,I_SW) = fact_ocean(i,j) * OCEAN_SFC_albedo(i,j,I_SW) &
                                         + fact_land (i,j) * LAND_SFC_albedo(i,j,I_SW) &
                                         + fact_urban(i,j) * URBAN_SFC_albedo(i,j,I_SW)
       ATMOS_PHY_SF_SFC_Z0M   (i,j)      = OCEAN_SFC_Z0M(i,j)
       ATMOS_PHY_SF_SFC_Z0H   (i,j)      = OCEAN_SFC_Z0H(i,j)
       ATMOS_PHY_SF_SFC_Z0E   (i,j)      = OCEAN_SFC_Z0E(i,j)
    end do
    end do


    !--- output boundary data
    totaltimesteps = totaltimesteps - NUMBER_OF_SKIP_TSTEPS ! skip first several data
    if ( totaltimesteps > 1 ) then
       if ( BOUNDARY_UPDATE_DT <= 0.0_RP ) then
          write(*,*) 'xxx BOUNDARY_UPDATE_DT is necessary in real case preprocess'
          call PRC_MPIstop
       endif

       call ParentOceanBoundary( OCEAN_TEMP_ORG(:,:,ns:ne),         &
                                 OCEAN_SFC_TEMP_ORG(:,:,ns:ne),     &
                                 OCEAN_SFC_albedo_ORG(:,:,:,ns:ne), &
                                 OCEAN_SFC_Z0_ORG(:,:,ns:ne),       &
                                 totaltimesteps,     &
                                 BOUNDARY_UPDATE_DT, &
                                 BASENAME_BOUNDARY,  &
                                 BOUNDARY_TITLE      )

    end if

    deallocate( ocean_temp_org )
    deallocate( ocean_sfc_temp_org )
    deallocate( ocean_sfc_albedo_org )
    deallocate( ocean_sfc_z0_org )

    return
  end subroutine REALINPUT_surface


  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtomSetup( &
      dims,                  &
      timelen,               &
      mdlid,                 &
      basename_org,          &
      filetype,              &
      use_file_density_in,   &
      serial_in              )
    use scale_external_io, only: &
         iSCALE, &
         iWRFARW, &
         iNICAM, &
         iGrADS
    use mod_realinput_scale, only: &
         ParentAtomSetupSCALE
    use mod_realinput_wrfarw, only: &
         ParentAtomSetupWRFARW
    use mod_realinput_nicam, only: &
         ParentAtomSetupNICAM
    use mod_realinput_grads, only: &
         ParentAtomSetupGrADS
    implicit none

    integer,          intent(out) :: dims(6)
    integer,          intent(out) :: timelen
    integer,          intent(out) :: mdlid
    character(len=*), intent(in)  :: basename_org
    character(len=*), intent(in)  :: filetype
    logical,          intent(in)  :: serial_in             ! read by a serial process
    logical,          intent(in)  :: use_file_density_in   ! use density data from files

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[RealinputAtmos]/Categ[Setup]'

    serial = serial_in
    if( serial ) then
       if( PRC_IsMaster ) then
          do_read_atom = .true.
       else
          do_read_atom = .false.
       endif
    else
       do_read_atom = .true.
    endif

    select case(trim(filetype))
    case('SCALE-RM')

       mdlid = iSCALE
       do_read_atom = .true.
       serial = .false.
       call ParentAtomSetupSCALE( dims ) ! (out)
       update_coord = .false.
       use_file_density = use_file_density_in
       use_temp = .false.
       rotate = .false.
       timelen = -1

    case('WRF-ARW')

       mdlid = iWRFARW
       if ( do_read_atom ) call ParentAtomSetupWRFARW( dims, timelen, & ! (out)
                                                       basename_org   ) ! (in)
       update_coord = .true.
       use_file_density = .false.
       use_temp = .true.
       rotate = .true.

    case('NICAM-NETCDF')

       mdlid = iNICAM
       if ( do_read_atom ) call ParentAtomSetupNICAM( dims, timelen, & ! (out)
                                                      basename_org   ) ! (in)
       update_coord = .false.
       use_file_density = .false.
       use_temp = .true.
       rotate = .true.

    case('GrADS')

       mdlid = iGrADS
       if ( do_read_atom ) call ParentAtomSetupGrADS( dims,        & ! (out)
                                                      basename_org ) ! (in)
       update_coord = .true.
       use_file_density = .false.
       use_temp = .true.
       rotate = .true.
       timelen = -1

    case default

       write(*,*) 'xxx Unsupported FILE TYPE:', trim(filetype)
       call PRC_MPIstop

    endselect

    if( serial ) then
       call COMM_bcast( dims(:), 6 )
       call COMM_bcast( timelen )
    endif

    if( IO_L ) write(IO_FID_LOG,*) '+++ Horizontal Interpolation Level:', &
                                    NEST_INTERP_LEVEL
    itp_nh = int( NEST_INTERP_LEVEL )
    itp_nv = 2

    allocate( hfact (        IA, JA, itp_nh         ) )
    allocate( vfact ( KA,    IA, JA, itp_nh, itp_nv ) )
    allocate( igrd  (        IA, JA, itp_nh         ) )
    allocate( jgrd  (        IA, JA, itp_nh         ) )
    allocate( kgrd  ( KA,    IA, JA, itp_nh, itp_nv ) )
    allocate( ncopy (        IA, JA, itp_nh         ) )

    allocate( lon_org (            dims(2), dims(3) ) )
    allocate( lat_org (            dims(2), dims(3) ) )
    allocate( cz_org  ( dims(1)+2, dims(2), dims(3) ) )

    allocate( velz_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( velx_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( vely_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( pott_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( temp_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( pres_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( qtrc_org( dims(1)+2, dims(2), dims(3), QA ) )
    allocate( dens_org( dims(1)+2, dims(2), dims(3) ) )

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
      dims,             &
      mdlid,            &
      flg_bin,          &
      flg_intrp,        &
      mptype_parent,    &
      timelen,          &
      skiplen           )
    use scale_comm, only: &
         COMM_vars8, &
         COMM_wait
    use scale_gridtrans, only: &
         rotc => GTRANS_ROTC
    use scale_atmos_thermodyn, only: &
         THERMODYN_pott => ATMOS_THERMODYN_pott
    use scale_atmos_hydrometeor, only: &
         HYDROMETEOR_diagnose_number_concentration => ATMOS_HYDROMETEOR_diagnose_number_concentration, &
         I_QV, &
         I_QC, &
         QLS, &
         QLE
    use scale_atmos_hydrostatic, only: &
         HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_interpolation_nest, only: &
         INTRPNEST_domain_compatibility, &
         INTRPNEST_interp_fact_llz, &
         INTRPNEST_interp_3d
    use mod_atmos_admin, only: &
         ATMOS_PHY_MP_TYPE
    use mod_realinput_scale, only: &
         ParentAtomOpenSCALE, &
         ParentAtomInputSCALE
    use mod_realinput_wrfarw, only: &
         ParentAtomOpenWRFARW, &
         ParentAtomInputWRFARW
    use mod_realinput_nicam, only: &
         ParentAtomOpenNICAM, &
         ParentAtomInputNICAM
    use mod_realinput_grads, only: &
         ParentAtomOpenGrADS, &
         ParentAtomInputGrADS
    implicit none

    real(RP),         intent(out) :: dens(:,:,:,:)
    real(RP),         intent(out) :: momz(:,:,:,:)
    real(RP),         intent(out) :: momx(:,:,:,:)
    real(RP),         intent(out) :: momy(:,:,:,:)
    real(RP),         intent(out) :: rhot(:,:,:,:)
    real(RP),         intent(out) :: qtrc(:,:,:,:,:)
    character(len=*), intent(in)  :: basename_org
    integer,          intent(in)  :: dims(6)
    integer,          intent(in)  :: mdlid            ! model type id
    logical,          intent(in)  :: flg_bin          ! flag for SBM(S10) is used or not 0-> not used, 1-> used
    logical,          intent(in)  :: flg_intrp        ! flag for interpolation of SBM(S10) from outer bulk-MP model
    integer,          intent(in)  :: mptype_parent    ! microphysics type of the parent model (number of classes)
    integer,          intent(in)  :: timelen          ! time steps in one file
    integer,          intent(in)  :: skiplen          ! skip steps

    real(RP) :: velz  (KA,IA,JA)
    real(RP) :: velx  (KA,IA,JA)
    real(RP) :: vely  (KA,IA,JA)
    real(RP) :: llvelx(KA,IA,JA)
    real(RP) :: llvely(KA,IA,JA)
    real(RP) :: work  (KA,IA,JA)
    real(RP) :: pott  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)

    real(RP) :: qc(KA,IA,JA)

    integer :: k, i, j, iq
    integer :: n, nn
    character(len=H_SHORT)  :: mptype_run
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[RealinputAtmos]/Categ[Input]'

    select case(ATMOS_PHY_MP_TYPE)
    case("DRY", "NONE")
       mptype_run = 'dry'
    case("KESSLER")
       mptype_run = 'single'
    case("TOMITA08")
       mptype_run = 'single'
    case("SN14")
       mptype_run = 'double'
    case("SUZUKI10")
       mptype_run = 'single-bin'
    case default
       write(*,*) 'xxx Unsupported ATMOS_PHY_MP_TYPE (', trim(ATMOS_PHY_MP_TYPE), '). Check!'
       call PRC_MPIstop
    end select


    if ( do_read_atom ) then

       select case( mdlid )
       case( iSCALE ) ! TYPE: SCALE-RM

          call ParentAtomOpenSCALE( lon_org, lat_org, & ! (out)
                                    cz_org,           & ! (out)
                                    basename_org,     & ! (in)
                                    dims              ) ! (in)

       case( iWRFARW ) ! TYPE: WRF-ARW

          call ParentAtomOpenWRFARW

       case( iNICAM ) ! TYPE: NICAM-NETCDF

          call ParentAtomOpenNICAM( lon_org, lat_org, & ! (out)
                                    cz_org,           & ! (out)
                                    basename_org,     & ! (in)
                                    dims              ) ! (in)

       case( iGrADS ) ! TYPE: GrADS format

          call ParentAtomOpenGrADS

       end select

    end if

    do n = skiplen+1, timelen
       nn = n - skiplen

       if ( do_read_atom ) then

          select case( mdlid )
          case( iSCALE ) ! TYPE: SCALE-RM

             call ParentAtomInputSCALE( velz_org, velx_org, vely_org, & ! (out)
                                        pres_org, dens_org, pott_org, & ! (out)
                                        qtrc_org,                     & ! (out)
                                        flg_bin, flg_intrp,           & ! (in)
                                        basename_org, mptype_parent,  & ! (in)
                                        dims, n         ) ! (in)

          case( iWRFARW ) ! TYPE: WRF-ARW

             call ParentAtomInputWRFARW( velz_org, velx_org, vely_org, & ! (out)
                                         pres_org, temp_org, qtrc_org, & ! (out)
                                         lon_org, lat_org, cz_org,     & ! (out)
                                         basename_org, mptype_parent,  & ! (in)
                                         dims, n                       ) ! (in)

          case( iNICAM ) ! TYPE: NICAM-NETCDF

             call ParentAtomInputNICAM( velz_org, velx_org, vely_org, & ! (out)
                                        pres_org, temp_org, qtrc_org, & ! (out)
                                        basename_org, dims, n         ) ! (in)

          case( iGrADS ) ! TYPE: GrADS format

             call ParentAtomInputGrADS( velz_org, velx_org, vely_org, & ! (out)
                                        pres_org, temp_org, qtrc_org, & ! (out)
                                        lon_org, lat_org, cz_org,     & ! (out)
                                        basename_org, dims, n         ) ! (in)

          end select

          if ( use_temp ) then
             do j = 1, dims(3)
             do i = 1, dims(2)
             do k = 1, dims(1)+2
                call THERMODYN_pott( pott_org(k,i,j),   & ! [OUT]
                                     temp_org(k,i,j),   & ! [IN]
                                     pres_org(k,i,j),   & ! [IN]
                                     qtrc_org(k,i,j,:), & ! [IN]
                                     TRACER_CV(:),      & ! [IN]
                                     TRACER_R(:),       & ! [IN]
                                     TRACER_MASS(:)     ) ! [IN]
             end do
             end do
             end do
          end if

       end if

       if ( serial ) then
          call COMM_bcast( velz_org, dims(1)+2, dims(2), dims(3) )
          call COMM_bcast( velx_org, dims(1)+2, dims(2), dims(3) )
          call COMM_bcast( vely_org, dims(1)+2, dims(2), dims(3) )
          call COMM_bcast( pott_org, dims(1)+2, dims(2), dims(3) )
          call COMM_bcast( qtrc_org, dims(1)+2, dims(2), dims(3), QA )
          if ( use_file_density ) then
             call COMM_bcast( dens_org, dims(1)+2, dims(2), dims(3) )
          else
             call COMM_bcast( pres_org, dims(1)+2, dims(2), dims(3) )
          end if

          if ( first .or. update_coord ) then

             call COMM_bcast( lon_org, dims(2), dims(3) )
             call COMM_bcast( lat_org, dims(2), dims(3) )
             call COMM_bcast( cz_org,  dims(1)+2, dims(2), dims(3) )

          end if

       end if

       if ( first .or. update_coord ) then

          call INTRPNEST_domain_compatibility( lon_org(:,:), lat_org(:,:), cz_org(:,:,:), &
                                               LON(:,:), LAT(:,:), CZ(KS:KE,:,:) )

          ! full level
          call INTRPNEST_interp_fact_llz( hfact, vfact,               & ! [OUT]
                                          kgrd, igrd, jgrd,           & ! [OUT]
                                          ncopy,                      & ! [OUT]
                                          CZ, LAT, LON,               & ! [IN]
                                          KS, KE, IA, JA,             & ! [IN]
                                          cz_org, lat_org, lon_org,   & ! [IN]
                                          dims(1)+2, dims(2), dims(3) ) ! [IN]

       end if

       call INTRPNEST_interp_3d( velz(:,:,:),      &
                                 velz_org(:,:,:),  &
                                 hfact(:,:,:),     &
                                 vfact(:,:,:,:,:), &
                                 kgrd(:,:,:,:,:),  &
                                 igrd(:,:,:),      &
                                 jgrd(:,:,:),      &
                                 IA, JA, KS, KE-1  )

       call INTRPNEST_interp_3d( llvelx  (:,:,:),     &
                                 velx_org(:,:,:),     &
                                 hfact   (:,:,:),     &
                                 vfact   (:,:,:,:,:), &
                                 kgrd    (:,:,:,:,:), &
                                 igrd    (:,:,:),     &
                                 jgrd    (:,:,:),     &
                                 IA, JA, KS, KE       )

       call INTRPNEST_interp_3d( llvely  (:,:,:),     &
                                 vely_org(:,:,:),     &
                                 hfact   (:,:,:),     &
                                 vfact   (:,:,:,:,:), &
                                 kgrd    (:,:,:,:,:), &
                                 igrd    (:,:,:),     &
                                 jgrd    (:,:,:),     &
                                 IA, JA, KS, KE       )

       if ( rotate ) then
          ! convert from latlon coordinate to local mapping (x)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             work(k,i,j) = llvelx(k,i,j) * rotc(i,j,cosin) + llvely(k,i,j) * rotc(i,j,sine )
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
          call COMM_wait ( velx(:,:,:), 1, .false. )

          ! convert from latlon coordinate to local mapping (y)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             work(k,i,j) = - llvelx(k,i,j) * rotc(i,j,sine ) + llvely(k,i,j) * rotc(i,j,cosin)
          end do
          end do
          end do

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
          call COMM_vars8( vely(:,:,:), 1 )
          call COMM_wait ( vely(:,:,:), 1, .false. )

       else
          velx = llvelx
          vely = llvely
       end if

       if( trim(mptype_run)=='double' .and. mptype_parent <= 6 )then
          if( IO_L ) write(IO_FID_LOG,*) '--- Diagnose Number Concentration from Mixing Ratio'
          call HYDROMETEOR_diagnose_number_concentration( qtrc_org(:,:,:,:) ) ! [inout]
       endif

       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)+2
          do iq = 1, QA
             qtrc_org(k,i,j,iq) = max( qtrc_org(k,i,j,iq), 0.0_RP )
          end do
       end do
       end do
       end do

       call INTRPNEST_interp_3d( pott    (:,:,:),       &
                                 pott_org(:,:,:),       &
                                 hfact   (:,:,:),     &
                                 vfact   (:,:,:,:,:), &
                                 kgrd    (:,:,:,:,:), &
                                 igrd    (:,:,:),     &
                                 jgrd    (:,:,:),     &
                                 IA, JA, KS, KE         )

       do iq = 1, QA
          call INTRPNEST_interp_3d( qtrc    (:,:,:,iq,nn), &
                                    qtrc_org(:,:,:,iq),    &
                                    hfact   (:,:,:),       &
                                    vfact   (:,:,:,:,:),   &
                                    kgrd    (:,:,:,:,:),   &
                                    igrd    (:,:,:),       &
                                    jgrd    (:,:,:),       &
                                    IA, JA, KS, KE         )
       end do

       if( use_file_density ) then
          ! use logarithmic density to interpolate more accurately

          dens_org = log( dens_org )
          call INTRPNEST_interp_3d( dens    (:,:,:,nn),  &
                                    dens_org(:,:,:),     &
                                    hfact   (:,:,:),     &
                                    vfact   (:,:,:,:,:), &
                                    kgrd    (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, KS, KE,      &
                                    logwegt=.true.       )
       else

          pres_org = log( pres_org )
          call INTRPNEST_interp_3d( pres    (:,:,:),     &
                                    pres_org(:,:,:),     &
                                    hfact   (:,:,:),     &
                                    vfact   (:,:,:,:,:), &
                                    kgrd    (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, KS, KE,      &
                                    logwegt=.true.       )

          qc = 0.0_RP
#ifndef DRY
          if ( I_QC > 0 ) then
             do iq = QLS, QLE
               qc(:,:,:) = qc(:,:,:) + QTRC(:,:,:,iq,nn)
             enddo
          end if
#endif
          ! make density & pressure profile in moist condition
          call HYDROSTATIC_buildrho_real( dens    (:,:,:,nn),      & ! [OUT]
                                          temp    (:,:,:),         & ! [OUT]
                                          pres    (:,:,:),         & ! [INOUT]
                                          pott    (:,:,:),         & ! [IN]
                                          qtrc    (:,:,:,I_QV,nn), & ! [IN]
                                          qc      (:,:,:)          ) ! [IN]

          call COMM_vars8( dens(:,:,:,nn), 1 )
          call COMM_wait ( dens(:,:,:,nn), 1 )

       end if


       do j = 1, JA
       do i = 1, IA
       do k = KS, KE-1
          momz(k,i,j,nn) = velz(k,i,j) * ( dens(k+1,i,j,nn) + dens(k,i,j,nn) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do i = 1, IA
          momz(KE,i,j,nn) = 0.0_RP
       end do
       end do
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          rhot(k,i,j,nn) = pott(k,i,j) * dens(k,i,j,nn)
       end do
       end do
       end do
       do j = 1, JA
       do i = 1, IA-1
       do k = KS, KE
          momx(k,i,j,nn) = velx(k,i,j) * ( dens(k,i+1,j,nn) + dens(k,i,j,nn) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do k = KS, KE
          momx(k,IA,j,nn) = velx(k,IA,j) * dens(k,IA,j,nn)
       end do
       end do
       call COMM_vars8( momx(:,:,:,nn), 1 )

       do j = 1, JA-1
       do i = 1, IA
       do k = KS, KE
          momy(k,i,j,nn) = vely(k,i,j) * ( dens(k,i,j+1,nn) + dens(k,i,j,nn) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          momy(k,i,JA,nn) = vely(k,i,JA) * dens(k,i,JA,nn)
       end do
       end do
       call COMM_vars8( momy(:,:,:,nn), 2 )

       call COMM_wait ( momx(:,:,:,nn), 1, .false. )
       call COMM_wait ( momy(:,:,:,nn), 2, .false. )

    end do

    first = .false.

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
       FILEIO_create, &
       FILEIO_def_var, &
       FILEIO_enddef, &
       FILEIO_write_var
    use scale_time, only: &
       TIME_NOWDATE
    use scale_atmos_phy_mp, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    implicit none

    real(RP),         intent(in)   :: dens(:,:,:,:)
    real(RP),         intent(in)   :: momz(:,:,:,:)
    real(RP),         intent(in)   :: momx(:,:,:,:)
    real(RP),         intent(in)   :: momy(:,:,:,:)
    real(RP),         intent(in)   :: rhot(:,:,:,:)
    real(RP),         intent(in)   :: qtrc(:,:,:,:,:)
    real(RP),         intent(in)   :: update_dt
    character(len=*), intent(in)   :: basename
    character(len=*), intent(in)   :: title
    integer,          intent(in)   :: numsteps ! total time steps

    character(len=H_SHORT) :: atmos_boundary_out_dtype = 'DEFAULT'  !< REAL4 or REAL8
    real(RP), allocatable :: buffer(:,:,:,:)
    integer :: nowdate(6)

    integer :: fid, vid(5+QA_MP)
    integer :: k, i, j, n, iq
    integer :: ts, te
    !---------------------------------------------------------------------------

    ts = 1
    te = numsteps

    nowdate = TIME_NOWDATE
    nowdate(1) = nowdate(1)

    allocate( buffer(KA,IA,JA,te-ts+1) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[RealinputAtmos]/Categ[Boundary]'

    call FILEIO_create( fid, basename, title, atmos_boundary_out_dtype, nowdate )

    call FILEIO_def_var( fid, vid(1), 'DENS', 'Reference Density', 'kg/m3', 'ZXYT', &
         atmos_boundary_out_dtype, update_dt, numsteps )
    call FILEIO_def_var( fid, vid(2), 'VELZ', 'Reference VELZ',    'm/s',   'ZXYT', &
         atmos_boundary_out_dtype, update_dt, numsteps )
    call FILEIO_def_var( fid, vid(3), 'VELX', 'Reference VELX',    'm/s',   'ZXYT', &
         atmos_boundary_out_dtype, update_dt, numsteps )
    call FILEIO_def_var( fid, vid(4), 'VELY', 'Reference VELY',    'm/s',   'ZXYT', &
         atmos_boundary_out_dtype, update_dt, numsteps )
    call FILEIO_def_var( fid, vid(5), 'POTT', 'Reference PT',      'K',     'ZXYT', &
         atmos_boundary_out_dtype, update_dt, numsteps )
    do iq = QS_MP, QE_MP
       call FILEIO_def_var( fid, vid(6+iq-QS_MP), TRACER_NAME(iq), 'Reference '//TRACER_NAME(iq), 'kg/kg', 'ZXYT', &
            atmos_boundary_out_dtype, update_dt, numsteps )
    end do

    call FILEIO_enddef( fid )

    call FILEIO_write_var( fid, vid(1), DENS(:,:,:,ts:te), 'DENS', 'ZXYT', update_dt )
    do n = ts, te
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMZ(k,i,j,n) / ( DENS(k+1,i,j,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    call FILEIO_write_var( fid, vid(2), buffer,            'VELZ', 'ZXYT', update_dt )
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
    call FILEIO_write_var( fid, vid(3), buffer,            'VELX', 'ZXYT', update_dt )
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
    call FILEIO_write_var( fid, vid(4), buffer,            'VELY', 'ZXYT', update_dt )
    do n = ts, te
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = RHOT(k,i,j,n) / DENS(k,i,j,n)
    end do
    end do
    end do
    end do
    call FILEIO_write_var( fid, vid(5), buffer,            'POTT', 'ZXYT', update_dt )

    do iq = QS_MP, QE_MP
       do n = ts, te
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          buffer(k,i,j,n-ts+1) = QTRC(k,i,j,iq,n)
       end do
       end do
       end do
       end do
       call FILEIO_write_var( fid, vid(6+iq-QS_MP), buffer, TRACER_NAME(iq), 'ZXYT', update_dt )
    end do

    deallocate( buffer )

    return
  end subroutine ParentAtomBoundary

  !-----------------------------------------------------------------------------
  !> Surface Setup
  subroutine ParentSurfaceSetup( &
       ldims, odims,        &
       lmdlid, omdlid,      &
       timelen,             &
       basename_land,       &
       basename_ocean,      &
       filetype_land,       &
       filetype_ocean,      &
       use_file_landwater,  &
       intrp_land_temp,     &
       intrp_land_water,    &
       intrp_land_sfc_temp, &
       intrp_ocean_temp,    &
       intrp_ocean_sfc_temp )
    use scale_external_io, only: &
         iSCALE, &
         iWRFARW, &
         iNICAM, &
         iGrADS
    use mod_realinput_scale, only: &
         ParentLandSetupSCALE, &
         ParentOceanSetupSCALE
    use mod_realinput_wrfarw, only: &
         ParentLandSetupWRFARW, &
         ParentOceanSetupWRFARW
    use mod_realinput_nicam, only: &
         ParentLandSetupNICAM, &
         ParentOceanSetupNICAM
    use mod_realinput_grads, only: &
         ParentLandSetupGrADS, &
         ParentOceanSetupGrADS
    implicit none

    integer,          intent(out) :: ldims(3) ! dims for land
    integer,          intent(out) :: odims(2) ! dims for ocean
    integer,          intent(out) :: lmdlid   ! model id for land
    integer,          intent(out) :: omdlid   ! model id for ocean
    integer,          intent(out) :: timelen  ! number of time steps in ocean file
    character(len=*), intent(in)  :: basename_land
    character(len=*), intent(in)  :: basename_ocean
    character(len=*), intent(in)  :: filetype_land
    character(len=*), intent(in)  :: filetype_ocean
    logical,          intent(in)  :: use_file_landwater ! use land water data from files
    character(len=*), intent(in)  :: intrp_land_temp
    character(len=*), intent(in)  :: intrp_land_water
    character(len=*), intent(in)  :: intrp_land_sfc_temp
    character(len=*), intent(in)  :: intrp_ocean_temp
    character(len=*), intent(in)  :: intrp_ocean_sfc_temp
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[RealinputSurface]/Categ[Setup]'

    ! Land

    if( LKMAX < 4 )then
       write(*,*) 'xxx LKMAX less than 4: ', LKMAX
       write(*,*) 'xxx in Real Case, LKMAX should be set more than 4'
       call PRC_MPIstop
    endif

    if( serial_land ) then
       if( PRC_IsMaster ) then
          do_read_land = .true.
       else
          do_read_land = .false.
       endif
    else
       do_read_land = .true.
    endif

    select case(trim(filetype_land))
    case('SCALE-RM')

       lmdlid = iSCALE
       serial_land = .false.
       do_read_land = .true.
       call ParentLandSetupSCALE( ldims ) ! (out)
       use_waterratio = .false.

    case('WRF-ARW')

       lmdlid = iWRFARW
       if ( do_read_land ) call ParentLandSetupWRFARW( ldims,        & ! (out)
                                                       basename_land ) ! (in)
       use_waterratio = .true.

    case('NICAM-NETCDF')

       lmdlid = iNICAM
       if ( do_read_land ) call ParentLandSetupNICAM( ldims,        & ! (out)
                                                      basename_land ) ! (in)
       use_waterratio = .false.

    case('GrADS')

       lmdlid = iGrADS
       if ( do_read_land ) call ParentLandSetupGrADS( ldims,              & ! (out)
                                                      use_waterratio,     & ! (out)
                                                      use_file_landwater, & ! (in)
                                                      basename_land       ) ! (in)

    case default

       write(*,*) 'xxx Unsupported FILE TYPE:', trim(filetype_land)
       call PRC_MPIstop

    endselect

    if( serial_land ) then
       call COMM_bcast( ldims(:), 3 )
       call COMM_bcast( use_waterratio )
    endif


    select case( INTRP_LAND_TEMP )
    case( 'off' )
       i_intrp_land_temp = i_intrp_off
    case( 'mask' )
       i_intrp_land_temp = i_intrp_mask
    case( 'fill' )
       i_intrp_land_temp = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_LAND_TEMP is invalid. ', INTRP_LAND_TEMP
       call PRC_MPIstop
    end select
    select case( INTRP_LAND_SFC_TEMP )
    case( 'off' )
       i_intrp_land_sfc_temp = i_intrp_off
    case( 'mask' )
       i_intrp_land_sfc_temp = i_intrp_mask
    case( 'fill' )
       i_intrp_land_sfc_temp = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_LAND_SFC_TEMP is invalid. ', INTRP_LAND_SFC_TEMP
       call PRC_MPIstop
    end select
    select case( INTRP_LAND_WATER )
    case( 'off' )
       i_intrp_land_water = i_intrp_off
    case( 'mask' )
       i_intrp_land_water = i_intrp_mask
    case( 'fill' )
       i_intrp_land_water = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_LAND_WATER is invalid. ', INTRP_LAND_WATER
       call PRC_MPIstop
    end select

    select case( lmdlid )
    case( iSCALE, iWRFARW, iNICAM )
       i_intrp_land_temp      = i_intrp_mask
       i_intrp_land_sfc_temp  = i_intrp_mask
       i_intrp_land_water     = i_intrp_mask
    end select


    ! Ocean

    if( serial_ocean ) then
       if( PRC_IsMaster ) then
          do_read_ocean = .true.
       else
          do_read_ocean = .false.
       endif
    else
       do_read_ocean = .true.
    endif

    select case(trim(filetype_ocean))
    case('SCALE-RM')

       timelen = -1
       omdlid = iSCALE
       serial_ocean = .false.
       do_read_ocean = .true.
       call ParentOceanSetupSCALE( odims )
       update_coord = .false.

    case('WRF-ARW')

       omdlid = iWRFARW
       if ( do_read_ocean ) call ParentOceanSetupWRFARW( odims, timelen, & ! (out)
                                                         basename_ocean  ) ! (in)
       update_coord = .true.

    case('NICAM-NETCDF')

       omdlid = iNICAM
       if ( do_read_ocean ) call ParentOceanSetupNICAM( odims, timelen, & ! (out)
                                                        basename_ocean  ) ! (in)
       update_coord = .false.

    case('GrADS')

       omdlid = iGrADS
       if ( do_read_ocean ) call ParentOceanSetupGrADS( odims, timelen, & ! (out)
                                                        basename_ocean  ) ! (out)
       update_coord = .false.

    case default

       write(*,*) 'xxx Unsupported FILE TYPE:', trim(filetype_ocean)
       call PRC_MPIstop

    endselect

    if( serial_ocean ) then
       call COMM_bcast( odims(:), 2 )
       call COMM_bcast( timelen )
    endif


    select case( INTRP_OCEAN_TEMP )
    case( 'off' )
       i_intrp_ocean_temp = i_intrp_off
    case( 'mask' )
       i_intrp_ocean_temp = i_intrp_mask
    case( 'fill' )
       i_intrp_ocean_temp = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_OCEAN_TEMP is invalid. ', INTRP_OCEAN_TEMP
       call PRC_MPIstop
    end select
    select case( INTRP_OCEAN_SFC_TEMP )
    case( 'off' )
       i_intrp_ocean_sfc_temp = i_intrp_off
    case( 'mask' )
       i_intrp_ocean_sfc_temp = i_intrp_mask
    case( 'fill' )
       i_intrp_ocean_sfc_temp = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_OCEAN_SFC_TEMP is invalid. ', INTRP_OCEAN_SFC_TEMP
       call PRC_MPIstop
    end select

    select case( omdlid )
    case( iSCALE, iWRFARW, iNICAM )
       i_intrp_ocean_temp     = i_intrp_mask
       i_intrp_ocean_sfc_temp = i_intrp_mask
    end select


    allocate( tw_org   ( odims(1), odims(2) ) )
    allocate( sst_org  ( odims(1), odims(2) ) )
    allocate( albw_org ( odims(1), odims(2), 2 ) )
    allocate( olon_org ( odims(1), odims(2) ) )
    allocate( olat_org ( odims(1), odims(2) ) )
    allocate( omask_org( odims(1), odims(2) ) )

    first = .true.

    return
  end subroutine ParentSurfaceSetup

  !-----------------------------------------------------------------------------
  !> Surface Data Read
  subroutine ParentSurfaceInput( &
       tg, &
       strg, &
       lst, &
       albg, &
       tc_urb, &
       qc_urb, &
       uc_urb, &
       ust, &
       albu, &
       tw, &
       sst, &
       albw, &
       z0w, &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       basename_land,     &
       basename_ocean,    &
       mdlid_land, &
       mdlid_ocean,            &
       ldims,             &
       odims,             &
       use_file_landwater, &
       init_landwater_ratio, &
       intrp_iter_max, &
       soilwater_ds2vc_flag, &
       elevation_collection, &
       timelen,          &
       skiplen, &
       lit )
    use scale_comm, only: &
         COMM_bcast, &
         COMM_vars8, &
         COMM_wait
    use scale_const, only: &
         EPS => CONST_EPS, &
         UNDEF => CONST_UNDEF, &
         I_SW => CONST_I_SW, &
         I_LW => CONST_I_LW
    use scale_interpolation_nest, only: &
         INTRPNEST_interp_fact_latlon, &
         INTRPNEST_interp_2d
    use scale_land_grid, only: &
         LCZ  => GRID_LCZ
    use scale_atmos_thermodyn, only: &
         THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_atmos_hydrometeor, only: &
         I_QV
    use scale_landuse, only: &
         lsmask_nest => LANDUSE_frac_land
    use mod_realinput_scale, only: &
         ParentOceanOpenSCALE, &
         ParentOCeanInputSCALE, &
         ParentLandInputSCALE
    use mod_realinput_wrfarw, only: &
         ParentOceanOpenWRFARW, &
         ParentOceanInputWRFARW, &
         ParentLandInputWRFARW
    use mod_realinput_nicam, only: &
         ParentOceanOpenNICAM, &
         ParentOceanInputNICAM, &
         ParentLandInputNICAM
    use mod_realinput_grads, only: &
         ParentOceanOpenGrADS, &
         ParentOceanInputGrADS, &
         ParentLandInputGrADS
    implicit none

    real(RP),         intent(inout) :: tg(LKMAX,IA,JA)
    real(RP),         intent(inout) :: strg(LKMAX,IA,JA)
    real(RP),         intent(inout) :: lst(IA,JA)
    real(RP),         intent(inout) :: albg(IA,JA,2)
    real(RP),         intent(inout) :: tc_urb(IA,JA)
    real(RP),         intent(inout) :: qc_urb(IA,JA)
    real(RP),         intent(inout) :: uc_urb(IA,JA)
    real(RP),         intent(inout) :: ust(IA,JA)
    real(RP),         intent(inout) :: albu(IA,JA,2)
    real(RP),         intent(out) :: tw(:,:,:)
    real(RP),         intent(out) :: sst(:,:,:)
    real(RP),         intent(out) :: albw(:,:,:,:)
    real(RP),         intent(out) :: z0w(:,:,:)
    real(RP),         intent(in)  :: DENS(KA,IA,JA)
    real(RP),         intent(in)  :: MOMZ(KA,IA,JA)
    real(RP),         intent(in)  :: MOMX(KA,IA,JA)
    real(RP),         intent(in)  :: MOMY(KA,IA,JA)
    real(RP),         intent(in)  :: RHOT(KA,IA,JA)
    real(RP),         intent(in)  :: QTRC(KA,IA,JA,QA)
    character(len=*), intent(in)  :: basename_land
    character(len=*), intent(in)  :: basename_ocean
    integer,          intent(in)  :: mdlid_land
    integer,          intent(in)  :: mdlid_ocean
    integer,          intent(in)  :: ldims(3)
    integer,          intent(in)  :: odims(2)
    logical,          intent(in)  :: use_file_landwater   ! use land water data from files
    real(RP),         intent(in)  :: init_landwater_ratio ! Ratio of land water to storage is constant,
                                                          ! if use_file_landwater is ".false."
    integer,          intent(in)  :: intrp_iter_max
    logical,          intent(in)  :: soilwater_ds2vc_flag
    logical,          intent(in)  :: elevation_collection
    integer,          intent(in)  :: timelen          ! time steps in one file
    integer,          intent(in)  :: skiplen          ! skip steps
    integer,          intent(in)  :: lit

   ! land
    real(RP) :: tg_org   (ldims(1),ldims(2),ldims(3))
    real(RP) :: strg_org (ldims(1),ldims(2),ldims(3))
    real(RP) :: smds_org (ldims(1),ldims(2),ldims(3))
!    real(RP) :: skint_org(         ldims(2),ldims(3))
    real(RP) :: lst_org  (         ldims(2),ldims(3))
    real(RP) :: ust_org  (         ldims(2),ldims(3))
    real(RP) :: albg_org (         ldims(2),ldims(3),2)
    real(RP) :: topo_org (         ldims(2),ldims(3))
    real(RP) :: lmask_org(         ldims(2),ldims(3))
    real(RP) :: lz_org   (ldims(1)                  )
    real(RP) :: llon_org (         ldims(2),ldims(3))
    real(RP) :: llat_org (         ldims(2),ldims(3))

    ! ocean
    real(RP) :: tw_org   (        odims(1),odims(2))
    real(RP) :: sst_org  (        odims(1),odims(2))
    real(RP) :: albw_org (        odims(1),odims(2),2)
    real(RP) :: z0w_org  (        odims(1),odims(2))
    real(RP) :: olon_org (        odims(1),odims(2))
    real(RP) :: olat_org (        odims(1),odims(2))
    real(RP) :: omask_org(        odims(1),odims(2))
    real(RP) :: omask    (        odims(1),odims(2))
    real(RP) :: lst_ocean(        odims(1),odims(2))

    real(RP) :: hfact_o(odims(1),odims(2),itp_nh)
    integer  :: igrd_o (odims(1),odims(2),itp_nh)
    integer  :: jgrd_o (odims(1),odims(2),itp_nh)

    real(RP) :: temp
    real(RP) :: pres

    integer :: i, j
    integer :: n, nn
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[RealinputOcean]/Categ[Input]'

    if ( first ) then ! read land data only once

       if ( do_read_land ) then

          select case( mdlid_land )
          case( iSCALE ) ! TYPE: SCALE-RM

             call ParentLandInputSCALE( &
                  tg_org, strg_org,           & ! (out)
                  lst_org, ust_org, albg_org, & ! (out)
                  topo_org, lmask_org,        & ! (out)
                  llon_org, llat_org, lz_org, & ! (out)
                  basename_land, ldims,       & ! (in)
                  use_file_landwater, lit     ) ! (in)

          case( iWRFARW ) ! TYPE: WRF-ARW

             call ParentLandInputWRFARW( &
                  tg_org, smds_org,           & ! (out)
                  lst_org, ust_org, albg_org, & ! (out)
                  topo_org, lmask_org,        & ! (out)
                  llon_org, llat_org, lz_org, & ! (out)
                  basename_land, ldims,       & ! (in)
                  use_file_landwater, lit     ) ! (in)

          case( iNICAM ) ! TYPE: NICAM-NETCDF

             call ParentLandInputNICAM( &
                  tg_org, strg_org,           & ! (out)
                  lst_org,                    & ! (out)
                  llon_org, llat_org, lz_org, & ! (out)
                  topo_org, lmask_org,        & ! (out)
                  basename_land, ldims,       & ! (in)
                  use_file_landwater, lit     ) ! (in)
             ust_org = UNDEF
             albg_org = UNDEF

          case( iGrADS ) ! TYPE: GrADS format

             call ParentLandInputGrADS( &
                  tg_org, strg_org, smds_org, & ! (out)
                  lst_org,                    & ! (out)
                  llon_org, llat_org, lz_org, & ! (out)
                  topo_org, lmask_org,        & ! (out)
                  basename_land, ldims,       & ! (in)
                  use_file_landwater, lit     ) ! (in)
             ust_org = UNDEF
             albg_org = UNDEF

          end select

       end if

       if ( serial_land ) then
          call COMM_bcast( tg_org, ldims(1), ldims(2), ldims(3) )
          if ( use_waterratio ) then
             call COMM_bcast( smds_org, ldims(1), ldims(2), ldims(3) )
          else
             call COMM_bcast( strg_org, ldims(1), ldims(2), ldims(3) )
          end if
          call COMM_bcast( lst_org, ldims(2), ldims(3) )
          call COMM_bcast( ust_org, ldims(2), ldims(3) )
          call COMM_bcast( albg_org(:,:,I_LW), ldims(2), ldims(3) )
          call COMM_bcast( albg_org(:,:,I_SW), ldims(2), ldims(3) )
          call COMM_bcast( topo_org, ldims(2), ldims(3) )
          call COMM_bcast( lmask_org, ldims(2), ldims(3) )
          call COMM_bcast( llon_org, ldims(2), ldims(3) )
          call COMM_bcast( llat_org, ldims(2), ldims(3) )
       end if


       ! urban data

       do j = 1, JA
       do i = 1, IA
          call THERMODYN_temp_pres( temp,           & ! [OUT]
                                    pres,           & ! [OUT] not used
                                    dens(KS,i,j),   & ! [IN]
                                    rhot(KS,i,j),   & ! [IN]
                                    qtrc(KS,i,j,:), & ! [IN]
                                    TRACER_CV(:),   & ! [IN]
                                    TRACER_R(:),    & ! [IN]
                                    TRACER_MASS(:)  ) ! [IN]

          tc_urb(i,j) = temp
#ifdef DRY
          qc_urb(i,j) = 0.0_RP
#else
          qc_urb(i,j) = qtrc(KS,i,j,I_QV)
#endif
       enddo
       enddo

       do j = 1, JA-1
       do i = 1, IA-1
          uc_urb(i,j) = max(sqrt( ( momx(KS,i,j) / (dens(KS,i+1,  j)+dens(KS,i,j)) * 2.0_RP )**2.0_RP &
                                + ( momy(KS,i,j) / (dens(KS,  i,j+1)+dens(KS,i,j)) * 2.0_RP )**2.0_RP ), &
                            0.01_RP)
       enddo
       enddo
       do j = 1, JA-1
          uc_urb(IA,j) = max(sqrt( ( momx(KS,IA,j) /  dens(KS,IA,j  ) )**2.0_RP &
                                 + ( momy(KS,IA,j) / (dens(KS,IA,j+1)+dens(KS,IA,j)) * 2.0_RP )**2.0_RP ), &
                             0.01_RP)
       enddo
       do i = 1, IA-1
          uc_urb(i,JA) = max(sqrt( ( momx(KS,i,JA) / (dens(KS,i+1,JA)+dens(KS,i,JA)) * 2.0_RP )**2.0_RP &
                                 + ( momy(KS,i,JA) /  dens(KS,i  ,JA) )**2.0_RP ), 0.01_RP)
       enddo
       uc_urb(IA,JA) = max(sqrt( ( momx(KS,IA,JA) / dens(KS,IA,JA) )**2.0_RP &
                               + ( momy(KS,IA,JA) / dens(KS,IA,JA) )**2.0_RP ), 0.01_RP)

       call COMM_vars8( uc_urb, 1 )
       call COMM_wait ( uc_urb, 1, .false. )


    end if ! first


    if ( do_read_ocean ) then

       select case( mdlid_ocean )
       case( iSCALE ) ! TYPE: SCALE-RM

          call ParentOceanOpenSCALE( olon_org, olat_org, & ! (out)
                                     omask_org,          & ! (out)
                                     basename_ocean,     & ! (in)
                                     odims               ) ! (in)

       case( iWRFARW ) ! TYPE: WRF-ARW

          call ParentOceanOpenWRFARW

       case( iNICAM ) ! TYPE: NICAM-NETCDF

          call ParentOceanOpenNICAM( olon_org, olat_org, & ! (out)
                                     omask_org,          & ! (out)
                                     basename_ocean,     & ! (in)
                                     odims               ) ! (in)

       case( iGrADS ) ! TYPE: GrADS format

          call ParentOceanOpenGrADS

       end select

    end if


    do n = skiplen+1, timelen
       nn = n - skiplen

       if ( do_read_ocean ) then

          select case( mdlid_ocean )
          case( iSCALE ) ! TYPE: SCALE-RM

             call ParentOceanInputSCALE( &
                  tw_org, sst_org,       & ! (out)
                  albw_org, z0w_org,     & ! (out)
                  omask_org,             & ! (out)
                  basename_ocean, odims, & ! (in)
                  n                      ) ! (in)

          case( iWRFARW ) ! TYPE: WRF-ARW

             call ParentOceanInputWRFARW( &
                  tw_org, sst_org,       & ! (out)
                  albw_org, z0w_org,     & ! (out)
                  omask_org,             & ! (out)
                  olon_org, olat_org,    & ! (out)
                  basename_ocean, odims, & ! (in)
                  n                      ) ! (in)

          case( iNICAM ) ! TYPE: NICAM-NETCDF

             call ParentOceanInputNICAM( &
                  tw_org, sst_org,       & ! (out)
                  basename_ocean, odims, & ! (in)
                  omask_org,             & ! (in)
                  n                      ) ! (in)
             albw_org = UNDEF
             z0w_org = UNDEF

          case( iGrADS ) ! TYPE: GrADS format

             call ParentOceanInputGrADS( &
                  tw_org, sst_org,       & ! (out)
                  omask_org,             & ! (out)
                  olon_org, olat_org,    & ! (out)
                  basename_ocean, odims, & ! (in)
                  n                      ) ! (in)
             albw_org = UNDEF
             z0w_org = UNDEF

          end select

       end if

       if ( serial_ocean ) then
          call COMM_bcast( tw_org, odims(1), odims(2) )
          call COMM_bcast( sst_org, odims(1), odims(2) )
          call COMM_bcast( albw_org(:,:,I_LW), odims(1), odims(2) )
          call COMM_bcast( albw_org(:,:,I_SW), odims(1), odims(2) )
          call COMM_bcast( z0w_org, odims(1), odims(2) )
          call COMM_bcast( omask_org, odims(1), odims(2) )
          if ( first .or. update_coord ) then
             call COMM_bcast( olon_org, odims(1), odims(2) )
             call COMM_bcast( olat_org, odims(1), odims(2) )
          end if
       end if


       if ( first .or. update_coord ) then

          ! interpolation facter between outer ocean grid
          call INTRPNEST_interp_fact_latlon( hfact_o(:,:,:),               & ! [OUT]
                                             igrd_o(:,:,:), jgrd_o(:,:,:), & ! [OUT]
                                             olat_org(:,:), olon_org(:,:), & ! [IN]
                                             odims(1), odims(2),           & ! [IN]
                                             llat_org(:,:), llon_org(:,:), & ! [IN]
                                             ldims(2), ldims(3)            ) ! [IN]

          call INTRPNEST_interp_fact_latlon( hfact_o(:,:,:),               & ! [OUT]
                                             igrd_o(:,:,:), jgrd_o(:,:,:), & ! [OUT]
                                             olat_org(:,:), olon_org(:,:), & ! [IN]
                                             odims(1), odims(2),           & ! [IN]
                                             llat_org(:,:), llon_org(:,:), & ! [IN]
                                             ldims(2), ldims(3)            ) ! [IN]

       end if

       ! Ocean temp: interpolate over the land
       if ( i_INTRP_OCEAN_TEMP .ne. i_intrp_off ) then
          select case( i_INTRP_OCEAN_TEMP )
          case( i_intrp_mask )
             omask = omask_org
          case( i_intrp_fill )
             call make_mask( omask, tw_org, odims(1), odims(2), landdata=.false.)
          end select
          call interp_OceanLand_data(tw_org, omask, odims(1), odims(2), .false., intrp_iter_max)
       end if

       ! SST: interpolate over the land
       if ( i_INTRP_OCEAN_SFC_TEMP .ne. i_intrp_off ) then
          select case( i_INTRP_OCEAN_SFC_TEMP )
          case( i_intrp_mask )
             omask = omask_org
          case( i_intrp_fill )
             call make_mask( omask, sst_org, odims(1), odims(2), landdata=.false.)
          end select
          call interp_OceanLand_data(sst_org, omask, odims(1), odims(2), .false., intrp_iter_max)
       end if

       if ( first ) then ! interporate land data only once

          call land_interporation( &
               tg, strg, & ! (out)
               lst, albg, & ! (out)
               ust, albu, & ! (out)
               tg_org, strg_org, smds_org, & ! (inout)
               lst_org, albg_org, & ! (inout)
               ust_org, & ! (inout)
               sst_org, & ! (in)
               lmask_org, & ! (in)
               lsmask_nest, & ! (in)
               topo_org, & ! (in)
               lz_org, llon_org, llat_org, & ! (in)
               LCZ, LON, LAT, & ! (in)
               ldims, odims, & ! (in)
               maskval_tg, maskval_strg, & ! (in)
               init_landwater_ratio, & ! (in)
               use_file_landwater, & ! (in)
               use_waterratio, & ! (in)
               soilwater_ds2vc_flag, & ! (in)
               elevation_collection, & ! (in)
               intrp_iter_max ) ! (in)

       end if ! first

       if ( first .or. update_coord ) then
          ! land surface temperature at ocean grid
          call INTRPNEST_interp_2d( lst_ocean(:,:), lst_org(:,:), hfact_o(:,:,:),        &
                                    igrd_o(:,:,:), jgrd_o(:,:,:), odims(1), odims(2) )
       end if

       call replace_misval_map( sst_org, lst_ocean, odims(1), odims(2), "SST")
       call replace_misval_map( tw_org,  lst_ocean, odims(1), odims(2), "OCEAN_TEMP")

       do j = 1, odims(2)
       do i = 1, odims(1)
          if ( albw_org(i,j,I_LW) == UNDEF ) albw_org(i,j,I_LW) = 0.04_RP  ! emissivity of water surface : 0.96
          if ( albw_org(i,j,I_SW) == UNDEF ) albw_org(i,j,I_SW) = 0.10_RP
          if ( z0w_org(i,j) == UNDEF ) z0w_org(i,j) = 0.001_RP
       end do
       end do


       if ( first .or. update_coord ) then
          ! interporation for ocean variables
          call INTRPNEST_interp_fact_latlon( hfact(:,:,:),                 & ! [OUT]
                                             igrd(:,:,:), jgrd(:,:,:),     & ! [OUT]
                                             LAT(:,:), LON(:,:),           & ! [IN]
                                             IA, JA,                       & ! [IN]
                                             olat_org(:,:), olon_org(:,:), & ! [IN]
                                             odims(1), odims(2)            ) ! [IN]
       end if

       call INTRPNEST_interp_2d( tw(:,:,nn), tw_org(:,:), hfact(:,:,:), &
                                 igrd(:,:,:), jgrd(:,:,:), IA, JA )
       call INTRPNEST_interp_2d( sst(:,:,nn), sst_org(:,:), hfact(:,:,:), &
                                 igrd(:,:,:), jgrd(:,:,:), IA, JA )
       call INTRPNEST_interp_2d( albw(:,:,I_LW,nn), albw_org(:,:,I_LW), hfact(:,:,:), &
                                 igrd(:,:,:), jgrd(:,:,:), IA, JA )
       call INTRPNEST_interp_2d( albw(:,:,I_SW,nn), albw_org(:,:,I_SW), hfact(:,:,:), &
                                 igrd(:,:,:), jgrd(:,:,:), IA, JA )
       call INTRPNEST_interp_2d( z0w(:,:,nn),   z0w_org(:,:),   hfact(:,:,:), &
                                 igrd(:,:,:), jgrd(:,:,:), IA, JA )

       if ( first ) then

          ! replace values over the ocean ####
          do j = 1, JA
          do i = 1, IA
             if( abs(lsmask_nest(i,j)-0.0_RP) < EPS ) then ! ocean grid
                lst(i,j)   = sst(i,j,nn)
                ust(i,j)   = sst(i,j,nn)
             endif
          enddo
          enddo

       end if


       first = .false.

    end do ! time loop

    return
  end subroutine ParentSurfaceInput

  !> Boundary Data Write
  subroutine ParentOceanBoundary( &
       tw, &
       sst, &
       albw, &
       z0, &
       numsteps,  &
       update_dt, &
       basename,  &
       title      )
    use scale_const, only: &
         I_SW => CONST_I_SW, &
         I_LW => CONST_I_LW
    use scale_fileio, only: &
         FILEIO_create, &
         FILEIO_def_var, &
         FILEIO_enddef, &
         FILEIO_write_var
    use scale_time, only: &
         TIME_NOWDATE
    implicit none

    real(RP),         intent(in)   :: tw(:,:,:)
    real(RP),         intent(in)   :: sst(:,:,:)
    real(RP),         intent(in)   :: albw(:,:,:,:)
    real(RP),         intent(in)   :: z0(:,:,:)
    real(RP),         intent(in)   :: update_dt
    character(len=*), intent(in)   :: basename
    character(len=*), intent(in)   :: title
    integer,          intent(in)   :: numsteps ! total time steps

    character(len=H_SHORT) :: ocean_boundary_out_dtype = 'DEFAULT'  !< REAL4 or REAL8
    integer :: nowdate(6)
    integer :: fid, vid(5)
    integer :: ts, te
    !---------------------------------------------------------------------------

    ts = 1
    te = numsteps

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[RealinputOcean]/Categ[Boundary]'

    nowdate = TIME_NOWDATE
    nowdate(1) = nowdate(1)

    call FILEIO_create( fid, basename, title, ocean_boundary_out_dtype, nowdate )

    call FILEIO_def_var( fid, vid(1), &
         'OCEAN_TEMP',     'Reference Ocean Temperature',            'K', 'XYT', &
         ocean_boundary_out_dtype, update_dt, numsteps )
    call FILEIO_def_var( fid, vid(2), &
         'OCEAN_SFC_TEMP', 'Reference Ocean Surface Temperature',    'K', 'XYT', &
         ocean_boundary_out_dtype, update_dt, numsteps )
    call FILEIO_def_var( fid, vid(3), &
         'OCEAN_ALB_LW', 'Reference Ocean Surface Albedo Long-wave', '1', 'XYT', &
         ocean_boundary_out_dtype, update_dt, numsteps )
    call FILEIO_def_var( fid, vid(4), &
         'OCEAN_ALB_SW', 'Reference Ocean Surface Albedo Short-wave', '1', 'XYT', &
         ocean_boundary_out_dtype, update_dt, numsteps )
    call FILEIO_def_var( fid, vid(5), &
         'OCEAN_SFC_Z0', 'Reference Ocean Surface Z0', 'm', 'XYT', &
         ocean_boundary_out_dtype, update_dt, numsteps )

    call FILEIO_enddef( fid )

       call FILEIO_write_var( fid, vid(1), tw(:,:,ts:te),         'OCEAN_TEMP',    'XYT', update_dt )
       call FILEIO_write_var( fid, vid(2), sst(:,:,ts:te),       'OCEAN_SFC_TEMP', 'XYT', update_dt )
       call FILEIO_write_var( fid, vid(3), albw(:,:,I_LW,ts:te), 'OCEAN_ALB_LW',   'XYT', update_dt )
       call FILEIO_write_var( fid, vid(4), albw(:,:,I_SW,ts:te), 'OCEAN_ALB_SW',   'XYT', update_dt )
       call FILEIO_write_var( fid, vid(5), z0(:,:,ts:te),        'OCEAN_SFC_Z0',   'XYT', update_dt )

    return
  end subroutine ParentOceanBoundary


  !-------------------------------
  subroutine land_interporation( &
       tg,                   &
       strg,                 &
       lst,                  &
       albg,                 &
       ust,                  &
       albu,                 &
       tg_org,               &
       strg_org,             &
       smds_org,             &
       lst_org,              &
       albg_org,             &
       ust_org,              &
       sst_org,              &
       lmask_org,            &
       lsmask_nest,          &
       topo_org,             &
       lz_org,               &
       llon_org,             &
       llat_org,             &
       LCZ,                  &
       LON,                  &
       LAT,                  &
       ldims,                &
       odims,                &
       maskval_tg,           &
       maskval_strg,         &
       init_landwater_ratio, &
       use_file_landwater,   &
       use_waterratio,       &
       soilwater_ds2vc_flag, &
       elevation_collection, &
       intrp_iter_max        )
    use scale_process, only: &
         PRC_MPIstop
    use scale_const, only: &
         UNDEF => CONST_UNDEF, &
         I_SW => CONST_I_SW, &
         I_LW => CONST_I_LW, &
         LAPS => CONST_LAPS
    use scale_interpolation_nest, only: &
         INTRPNEST_interp_fact_llz, &
         INTRPNEST_interp_fact_latlon, &
         INTRPNEST_interp_3d, &
         INTRPNEST_interp_2d
    use scale_topography, only: &
         TOPO_Zsfc
    use mod_land_vars, only: &
         convert_WS2VWC
    implicit none
    real(RP), intent(out)   :: tg(LKMAX,IA,JA)
    real(RP), intent(out)   :: strg(LKMAX,IA,JA)
    real(RP), intent(out)   :: lst(IA,JA)
    real(RP), intent(out)   :: albg(IA,JA,2)
    real(RP), intent(out)   :: ust(IA,JA)
    real(RP), intent(out)   :: albu(IA,JA,2)
    real(RP), intent(inout) :: tg_org(:,:,:)
    real(RP), intent(inout) :: strg_org(:,:,:)
    real(RP), intent(inout) :: smds_org(:,:,:)
    real(RP), intent(inout) :: lst_org(:,:)
    real(RP), intent(inout) :: albg_org(:,:,:)
    real(RP), intent(inout) :: ust_org(:,:)
    real(RP), intent(inout) :: sst_org(:,:)
    real(RP), intent(in)    :: lmask_org(:,:)
    real(RP), intent(in)    :: lsmask_nest(:,:)
    real(RP), intent(in)    :: topo_org(:,:)
    real(RP), intent(in)    :: lz_org(:)
    real(RP), intent(in)    :: llon_org(:,:)
    real(RP), intent(in)    :: llat_org(:,:)
    real(RP), intent(in)    :: LCZ(LKMAX)
    real(RP), intent(in)    :: LON(IA,JA)
    real(RP), intent(in)    :: LAT(IA,JA)
    integer,  intent(in)    :: ldims(3)
    integer,  intent(in)    :: odims(2)
    real(RP), intent(in)    :: maskval_tg
    real(RP), intent(in)    :: maskval_strg
    real(RP), intent(in)    :: init_landwater_ratio
    logical,  intent(in)    :: use_file_landwater
    logical,  intent(in)    :: use_waterratio
    logical,  intent(in)    :: soilwater_ds2vc_flag
    logical,  intent(in)    :: elevation_collection
    integer,  intent(in)    :: intrp_iter_max

    real(RP) :: lmask(ldims(2), ldims(3))
    real(RP) :: smds(LKMAX,IA,JA)

    ! data for interporation
    real(RP) :: hfact_l(ldims(2), ldims(3), itp_nh)
    integer  :: igrd_l (ldims(2), ldims(3), itp_nh)
    integer  :: jgrd_l (ldims(2), ldims(3), itp_nh)
    real(RP) :: vfactl(LKMAX,IA,JA,itp_nh,itp_nv)
    integer  :: kgrdl (LKMAX,IA,JA,itp_nh,itp_nv)

    real(RP) :: sst_land(ldims(2), ldims(3))
    real(RP) :: work(ldims(2), ldims(3))

    real(RP) :: lz3d_org(ldims(1),ldims(2),ldims(3))
    real(RP) :: lcz_3D(LKMAX,IA,JA)

    ! elevation collection
    real(RP) :: topo(IA,JA)
    real(RP) :: tdiff

    integer :: k, i, j


    ! Surface skin temp: interpolate over the ocean
    if ( i_INTRP_LAND_SFC_TEMP .ne. i_intrp_off ) then
       select case( i_INTRP_LAND_SFC_TEMP )
       case( i_intrp_mask )
          lmask = lmask_org
       case( i_intrp_fill )
          call make_mask( lmask, lst_org, ldims(2), ldims(3), landdata=.true.)
       case default
          write(*,*) 'xxx INTRP_LAND_SFC_TEMP is invalid.'
          call PRC_MPIstop
       end select
       call interp_OceanLand_data(lst_org, lmask, ldims(2), ldims(3), .true., intrp_iter_max)
    end if

    ! Urban surface temp: interpolate over the ocean
    ! if ( i_INTRP_URB_SFC_TEMP .ne. i_intrp_off ) then
    !   select case( i_INTRP_URB_SFC_TEMP )
    !   case( i_intrp_mask )
    !      lmask = lmask_org
    !   case( i_intrp_fill )
    !      call make_mask( lmask, ust_org, ldims(2), ldims(3), landdata=.true.)
    !   case default
    !      write(*,*) 'xxx INTRP_URB_SFC_TEMP is invalid.'
    !      call PRC_MPIstop
    !   end select
    !   call interp_OceanLand_data(ust_org, lmask, ldims(2), ldims(3), .true., intrp_iter_max)
    !end if

    ! interpolation facter between outer land grid and ocean grid
    call INTRPNEST_interp_fact_latlon( hfact_l(:,:,:),               & ! [OUT]
                                       igrd_l(:,:,:), jgrd_l(:,:,:), & ! [OUT]
                                       llat_org(:,:), llon_org(:,:), & ! [IN]
                                       ldims(2), ldims(3),           & ! [IN]
                                       olat_org(:,:), olon_org(:,:), & ! [IN]
                                       odims(1), odims(2)            ) ! [IN]


    ! sst on land grid
    call INTRPNEST_interp_2d( sst_land(:,:), sst_org(:,:), hfact_l(:,:,:),      &
                              igrd_l(:,:,:), jgrd_l(:,:,:), ldims(2), ldims(3) )
    call replace_misval_map( lst_org, sst_land, ldims(2), ldims(3), "SKINT")

    ! replace missing value
    do j = 1, ldims(3)
    do i = 1, ldims(2)
       if ( ust_org(i,j) == UNDEF ) ust_org(i,j) = lst_org(i,j)
!       if ( skinw_org(i,j) == UNDEF ) skinw_org(i,j) = 0.0_RP
!       if ( snowq_org(i,j) == UNDEF ) snowq_org(i,j) = 0.0_RP
!       if ( snowt_org(i,j) == UNDEF ) snowt_org(i,j) = TEM00
       if ( albg_org(i,j,I_LW) == UNDEF ) albg_org(i,j,I_LW) = 0.03_RP  ! emissivity of general ground surface : 0.95-0.98
       if ( albg_org(i,j,I_SW) == UNDEF ) albg_org(i,j,I_SW) = 0.22_RP
    end do
    end do

    ! Land temp: interpolate over the ocean
    if ( i_INTRP_LAND_TEMP .ne. i_intrp_off ) then
       do k = 1, ldims(1)
          work(:,:) = tg_org(k,:,:)
          select case( i_INTRP_LAND_TEMP )
          case( i_intrp_mask )
             lmask = lmask_org
          case( i_intrp_fill )
             call make_mask( lmask, work, ldims(2), ldims(3), landdata=.true.)
          end select
          call interp_OceanLand_data( work, lmask, ldims(2), ldims(3), .true., intrp_iter_max )
          !replace land temp using skin temp
          call replace_misval_map( work, lst_org, ldims(2),  ldims(3),  "STEMP")
          tg_org(k,:,:) = work(:,:)
       end do
    end if


    ! fill grid data
    do j = 1, ldims(3)
    do i = 1, ldims(2)
       lz3d_org(:,i,j) = lz_org(:)
    end do
    end do

    do j = 1, JA
    do i = 1, IA
       lcz_3D(:,i,j) = LCZ(:)
    enddo
    enddo

    call INTRPNEST_interp_fact_llz( hfact  (:,:,:),               & ! [OUT]
                                    vfactl (:,:,:,:,:),           & ! [OUT]
                                    kgrdl  (:,:,:,:,:),           & ! [OUT]
                                    igrd   (:,:,:), jgrd(:,:,:),  & ! [OUT]
                                    ncopy  (:,:,:),               & ! [OUT]
                                    lcz_3D (:,:,:),               & ! [IN]
                                    LAT    (:,:), LON    (:,:),   & ! [IN]
                                    1, LKMAX, IA, JA,             & ! [IN]
                                    lz3d_org(:,:,:),              & ! [IN]
                                    llat_org(:,:), llon_org(:,:), & ! [IN]
                                    ldims(1), ldims(2), ldims(3), & ! [IN]
                                    landgrid=.true.               ) ! [IN]

    call INTRPNEST_interp_2d( lst(:,:),   lst_org(:,:),   hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
    call INTRPNEST_interp_2d( ust(:,:),   ust_org(:,:),   hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
!    call INTRPNEST_interp_2d( skinw(:,:), skinw_org(:,:), hfact(:,:,:), &
!                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
!    call INTRPNEST_interp_2d( snowq(:,:), snowq_org(:,:), hfact(:,:,:), &
!                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
!    call INTRPNEST_interp_2d( snowt(:,:), snowt_org(:,:), hfact(:,:,:), &
!                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
    call INTRPNEST_interp_2d( albg(:,:,I_LW), albg_org(:,:,I_LW), hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
    call INTRPNEST_interp_2d( albg(:,:,I_SW), albg_org(:,:,I_SW), hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )

    call INTRPNEST_interp_3d( tg    (:,:,:),     &
                              tg_org(:,:,:),     &
                              hfact (:,:,:),     &
                              vfactl(:,:,:,:,:), &
                              kgrdl (:,:,:,:,:), &
                              igrd  (:,:,:),     &
                              jgrd  (:,:,:),     &
                              IA, JA, 1, LKMAX-1 )

    do j = 1, JA
    do i = 1, IA
       tg(LKMAX,i,j) = tg(LKMAX-1,i,j)
    enddo ! i
    enddo ! j

    ! replace values over the ocean
    do k = 1, LKMAX
       call replace_misval_const( tg(k,:,:), maskval_tg, lsmask_nest )
    enddo


    ! elevation collection
    if ( elevation_collection ) then
       call INTRPNEST_interp_2d( topo(:,:),   topo_org(:,:),   hfact(:,:,:), &
                                 igrd(:,:,:), jgrd(:,:,:), IA, JA )

       do j = 1, JA
       do i = 1, IA
          if ( topo(i,j) > 0.0_RP ) then ! ignore UNDEF value
             tdiff = ( TOPO_Zsfc(i,j) - topo(i,j) ) * LAPS
             lst(i,j) = lst(i,j) - tdiff
             ust(i,j) = ust(i,j) - tdiff
             do k = 1, LKMAX
                tg(k,i,j) = tg(k,i,j) - tdiff
             end do
          end if
       end do
       end do
    end if



    ! Land water: interpolate over the ocean
    if( use_file_landwater )then

       if ( use_waterratio ) then

          if ( i_INTRP_LAND_WATER .ne. i_intrp_off ) then
             do k = 1, ldims(1)
                work(:,:) = smds_org(k,:,:)
                select case( i_INTRP_LAND_WATER )
                case( i_intrp_mask )
                   lmask = lmask_org
                case( i_intrp_fill )
                   call make_mask( lmask, work, ldims(2), ldims(3), landdata=.true.)
                end select
                call interp_OceanLand_data(work, lmask, ldims(2), ldims(3), .true., intrp_iter_max)
                lmask(:,:) = init_landwater_ratio
                !replace missing value to init_landwater_ratio
                call replace_misval_map( work, lmask, ldims(2), ldims(3),  "SMOISDS")
                smds_org(k,:,:) = work(:,:)
             enddo
          end if

          call INTRPNEST_interp_3d( smds    (:,:,:),     &
                                    smds_org(:,:,:),     &
                                    hfact   (:,:,:),     &
                                    vfactl  (:,:,:,:,:), &
                                    kgrdl   (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, 1, LKMAX-1   )
          do k = 1, LKMAX
             strg(k,:,:) = convert_WS2VWC( smds(k,:,:), critical=soilwater_DS2VC_flag )
          end do

       else

          if ( i_INTRP_LAND_WATER .ne. i_intrp_off ) then
             do k = 1, ldims(1)
                work(:,:) = strg_org(k,:,:)
                select case( i_INTRP_LAND_WATER )
                case( i_intrp_mask )
                   lmask = lmask_org
                case( i_intrp_fill )
                   call make_mask( lmask, work, ldims(2), ldims(3), landdata=.true.)
                end select
                call interp_OceanLand_data(work, lmask, ldims(2), ldims(3), .true., intrp_iter_max)
                lmask(:,:) = maskval_strg
                !replace missing value to init_landwater_ratio
                call replace_misval_map( work, lmask, ldims(2), ldims(3),  "SMOIS")
                strg_org(k,:,:) = work(:,:)
             enddo
          end if

          call INTRPNEST_interp_3d( strg    (:,:,:),     &
                                    strg_org(:,:,:),     &
                                    hfact   (:,:,:),     &
                                    vfactl  (:,:,:,:,:), &
                                    kgrdl   (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, 1, LKMAX-1   )
          ! interpolation
          do j = 1, JA
          do i = 1, IA
             strg(LKMAX,i,j) = strg(LKMAX-1,i,j)
          enddo
          enddo

       end if

       ! replace values over the ocean
       do k = 1, LKMAX
          call replace_misval_const( strg(k,:,:), maskval_strg, lsmask_nest )
       enddo

    else  ! not read from boundary file

       smds(:,:,:) = init_landwater_ratio
       ! conversion from water saturation [fraction] to volumetric water content [m3/m3]
       do k = 1, LKMAX
          strg(k,:,:) = convert_WS2VWC( smds(k,:,:), critical=.true. )
       end do

    endif ! use_file_waterratio


    ! copy albedo of land to urban
    do j = 1, JA
    do i = 1, IA
       albu(i,j,:) = albg(i,j,:)
    enddo
    enddo


    return
  end subroutine land_interporation

  !-------------------------------
  subroutine make_mask( &
      gmask,     & ! (out)
      data,      & ! (in)
      nx,        & ! (in)
      ny,        & ! (in)
      landdata   ) ! (in)
    use scale_const, only: &
       EPS => CONST_EPS, &
       UNDEF => CONST_UNDEF
    implicit none
    real(RP), intent(out)  :: gmask(:,:)
    real(RP), intent(in)   :: data(:,:)
    integer,  intent(in)   :: nx
    integer,  intent(in)   :: ny
    logical,  intent(in)   :: landdata   ! .true. => land data , .false. => ocean data

    real(RP)               :: dd
    integer                :: i,j

    if( landdata )then
       gmask(:,:) = 1.0_RP  ! gmask=1 will be skip in "interp_OceanLand_data"
       dd         = 0.0_RP
    else
       gmask(:,:) = 0.0_RP  ! gmask=0 will be skip in "interp_OceanLand_data"
       dd         = 1.0_RP
    endif

    do j = 1, ny
    do i = 1, nx
       if( abs(data(i,j) - UNDEF) < sqrt(EPS) )then
          gmask(i,j) = dd
       endif
    enddo
    enddo

    return
  end subroutine make_mask
  !-----------------------------------------------------------------------------
  subroutine interp_OceanLand_data( &
      data,      & ! (inout)
      lsmask,    & ! (in)
      nx,        & ! (in)
      ny,        & ! (in)
      landdata,  & ! (in)
      iter_max,  & ! (in)
      maskval    & ! (out)
      )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none
    real(RP), intent(inout)         :: data(:,:)
    real(RP), intent(in)            :: lsmask(:,:)
    integer,  intent(in)            :: nx
    integer,  intent(in)            :: ny
    logical,  intent(in)            :: landdata   ! .true. => land data , .false. => ocean data
    integer,  intent(in)            :: iter_max
    real(RP), intent(out), optional :: maskval

    integer                 :: untarget_mask
    integer, allocatable    :: imask(:,:),imaskr(:,:)
    real(RP),allocatable    :: newdata(:,:)
    real(RP)                :: nd
    integer                 :: count

    integer :: i, j, ii, jj, kk

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

    imaskr = imask
    do kk = 1, iter_max
       do j  = 1, ny
       do i  = 1, nx
          if ( imask(i,j) == untarget_mask ) then  ! not missing value
             newdata(i,j) = data(i,j)
             cycle
          else

             if ( present(maskval) ) then
             if ( abs(maskval-999.99_RP)<EPS ) then
                if (abs(lsmask(i,j)-0.0_RP) < EPS) maskval = data(i,j)
             endif
             endif

             !--------------------------------------
             ! check data of neighbor grid
             !---------------------------------------
             count = 0
             nd = 0.0_RP
             do  jj = j-1, j+1
                if ( jj < 1 .or. jj > ny ) cycle
                do ii = i-1, i+1
                   if ( ii < 1 .or. ii > nx .or. (jj == j .and. ii == i) ) cycle
                   if ( imask(ii,jj) == untarget_mask ) then
                      nd = nd + data(ii,jj)
                      count = count + 1
                   end if
                end do
             end do

             if( count >= 3 )then  ! coast grid : interpolate
                newdata(i,j) = nd / count
                imaskr(i,j) = untarget_mask
             else
                newdata(i,j) = data(i,j)
             endif

          endif ! sea/land

       enddo
       enddo

       imask(:,:) = imaskr(:,:)
       data(:,:)  = newdata(:,:)
    enddo ! kk

    deallocate( imask   )
    deallocate( imaskr  )
    deallocate( newdata )

    return
  end subroutine interp_OceanLand_data

  !-----------------------------------------------------------------------------
  subroutine replace_misval_const( data, maskval, frac_land )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none
    real(RP), intent(inout) :: data(:,:)
    real(RP), intent(in)    :: maskval
    real(RP), intent(in)    :: frac_land(:,:)
    integer                 :: i, j

    do j = 1, JA
    do i = 1, IA
       if( abs(frac_land(i,j)-0.0_RP) < EPS )then ! ocean grid
          data(i,j) = maskval
       endif
    enddo
    enddo

  end subroutine replace_misval_const

  !-----------------------------------------------------------------------------
  subroutine replace_misval_map( data, maskval, nx, ny, elem)
    use scale_const, only: &
       EPS => CONST_EPS, &
       UNDEF => CONST_UNDEF
    implicit none

    real(RP),         intent(inout) :: data(:,:)
    real(RP),         intent(in)    :: maskval(:,:)
    integer,          intent(in)    :: nx, ny
    character(len=*), intent(in)    :: elem

    integer :: i, j

    do j = 1, ny
    do i = 1, nx
       if( abs(data(i,j) - UNDEF) < sqrt(EPS) )then
          if( abs(maskval(i,j) - UNDEF) < sqrt(EPS) )then
             write(*,*) "Data for mask has missing value. ",trim(elem),i,j
             call PRC_MPIstop
          else
             data(i,j) = maskval(i,j)
          endif
       endif
    enddo
    enddo

  end subroutine replace_misval_map

end module mod_realinput
