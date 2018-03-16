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
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_ocean_grid_cartesC_index
  use scale_land_grid_cartesC_index
  use scale_urban_grid_cartesC_index
  use scale_index
  use scale_tracer

  use scale_process, only: &
     PRC_IsMaster, &
     PRC_MPIstop
  use scale_comm, only: &
     COMM_bcast
  use scale_atmos_grid_cartesC_real, only: &
     LON => ATMOS_GRID_CARTESC_REAL_LON, &
     LAT => ATMOS_GRID_CARTESC_REAL_LAT, &
     CZ  => ATMOS_GRID_CARTESC_REAL_CZ,  &
     FZ  => ATMOS_GRID_CARTESC_REAL_FZ
  use scale_comm_cartesC_nest, only: &
     COMM_CARTESC_NEST_INTERP_LEVEL
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
  private :: ParentAtmosSetup
  private :: ParentAtmosOpen
  private :: ParentAtmosInput
  private :: BoundaryAtmosSetup
  private :: BoundaryAtmosOutput

  private :: ParentSurfaceSetup
  private :: ParentSurfaceInput
  private :: ParentSurfaceBoundary
  private :: interp_OceanLand_data

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, public, parameter :: iSCALE  = 1
  integer, public, parameter :: iWRFARW = 2
!  integer, public, parameter :: iNICAM  = 3
  integer, public, parameter :: iGrADS  = 4

  real(RP), private, allocatable :: LON_org (:,:)
  real(RP), private, allocatable :: LAT_org (:,:)
  real(RP), private, allocatable :: CZ_org  (:,:,:)

  real(RP), private, allocatable :: W_org   (:,:,:) ! scalar point
  real(RP), private, allocatable :: U_org   (:,:,:) ! scalar point
  real(RP), private, allocatable :: V_org   (:,:,:) ! scalar point
  real(RP), private, allocatable :: DENS_org(:,:,:)
  real(RP), private, allocatable :: POTT_org(:,:,:)
  real(RP), private, allocatable :: TEMP_org(:,:,:)
  real(RP), private, allocatable :: PRES_org(:,:,:)
  real(RP), private, allocatable :: QTRC_org(:,:,:,:)
  real(RP), private, allocatable :: QV_org  (:,:,:)
  real(RP), private, allocatable :: QHYD_org(:,:,:,:)
  real(RP), private, allocatable :: QNUM_org(:,:,:,:)

  integer,  private, allocatable :: igrd (:,:,:)
  integer,  private, allocatable :: jgrd (:,:,:)
  real(RP), private, allocatable :: hfact(:,:,:)
  integer,  private, allocatable :: kgrd (:,:,:,:,:)
  real(RP), private, allocatable :: vfact(:,:,:,:,:)

  real(RP), private, allocatable :: tw_org   (:,:)
  real(RP), private, allocatable :: sst_org  (:,:)
  real(RP), private, allocatable :: albw_org (:,:,:)
  real(RP), private, allocatable :: olon_org (:,:)
  real(RP), private, allocatable :: olat_org (:,:)
  real(RP), private, allocatable :: omask_org(:,:)

  integer,  private              :: itp_nh = 4
  integer,  private              :: itp_nv = 2

  logical,  private              :: serial_atmos
  logical,  private              :: serial_land
  logical,  private              :: serial_ocean
  logical,  private              :: read_by_myproc_atmos
  logical,  private              :: do_read_land
  logical,  private              :: do_read_ocean

  logical,  private              :: temp2pott
  logical,  private              :: apply_rotate_uv
  logical,  private              :: update_coord
  logical,  private              :: use_waterratio

  integer,  private, parameter   :: I_intrp_off  = 0
  integer,  private, parameter   :: I_intrp_mask = 1
  integer,  private, parameter   :: I_intrp_fill = 2

  integer,  private              :: i_intrp_land_temp
  integer,  private              :: i_intrp_land_water
  integer,  private              :: i_intrp_land_sfc_temp
  integer,  private              :: i_intrp_ocean_temp
  integer,  private              :: i_intrp_ocean_sfc_temp

  ! replace missing value
  real(RP), private, parameter   :: maskval_tg   = 298.0_RP ! mask value 298K
  real(RP), private, parameter   :: maskval_strg = 0.02_RP  ! mask value 0.02
                                                            ! default value 0.02: set as value of forest at 40% of evapolation rate.
                                                            ! forest is considered as a typical landuse over Japan area.

  ! for namelist
  integer,                private :: NUMBER_OF_FILES            = 1
  integer,                private :: NUMBER_OF_TSTEPS           = 1       ! num of time steps in one file
  integer,                private :: NUMBER_OF_SKIP_TSTEPS      = 0       ! num of skipped first several data

  logical,                private :: SERIAL_PROC_READ           = .true.  ! read by one MPI process and broadcast

  character(len=H_LONG),  private :: FILETYPE_ORG               = ''
  character(len=H_LONG),  private :: BASENAME_ORG               = ''
  logical,                private :: BASENAME_ADD_NUM           = .false.

  character(len=H_LONG),  private :: BASENAME_BOUNDARY          = ''
  logical,                private :: BOUNDARY_POSTFIX_TIMELABEL = .false.
  character(len=H_LONG),  private :: BOUNDARY_TITLE             = 'SCALE-RM BOUNDARY CONDITION for REAL CASE'
  character(len=H_SHORT), private :: BOUNDARY_DTYPE             = 'DEFAULT'
  real(DP),               private :: BOUNDARY_UPDATE_DT         = 0.0_DP  ! inteval time of boudary data update [s]

  logical,                private :: USE_FILE_DENSITY           = .false. ! use density data from files
  logical,                private :: SAME_MP_TYPE               = .false. ! microphysics type of the parent model is same as it in this model

  logical, private :: first = .true.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine REALINPUT_atmos
    use scale_time, only: &
       TIME_gettimelabel
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

    NAMELIST / PARAM_MKINIT_REAL_ATMOS / &
       NUMBER_OF_FILES,            &
       NUMBER_OF_TSTEPS,           &
       NUMBER_OF_SKIP_TSTEPS,      &
       SERIAL_PROC_READ,           &
       FILETYPE_ORG,               &
       BASENAME_ORG,               &
       BASENAME_ADD_NUM,           &
       BASENAME_BOUNDARY,          &
       BOUNDARY_POSTFIX_TIMELABEL, &
       BOUNDARY_TITLE,             &
       BOUNDARY_DTYPE,             &
       BOUNDARY_UPDATE_DT,         &
       USE_FILE_DENSITY,           &
       SAME_MP_TYPE

    character(len=H_LONG) :: basename_mod
    character(len=H_LONG) :: basename_out_mod
    character(len=19)     :: timelabel

    integer  :: dims(6) ! dims 1-3: normal, 4-6: staggerd
    integer  :: timelen

    integer  :: fid_atmos
    integer  :: vid_atmos(5+QA)

    real(RP) :: DENS_in(KA,IA,JA)
    real(RP) :: MOMZ_in(KA,IA,JA) ! staggered point
    real(RP) :: MOMX_in(KA,IA,JA) ! staggered point
    real(RP) :: MOMY_in(KA,IA,JA) ! staggered point
    real(RP) :: RHOT_in(KA,IA,JA)
    real(RP) :: QTRC_in(KA,IA,JA,QA)

    real(RP) :: VELZ_in(KA,IA,JA) ! staggered point
    real(RP) :: VELX_in(KA,IA,JA) ! staggered point
    real(RP) :: VELY_in(KA,IA,JA) ! staggered point
    real(RP) :: POTT_in(KA,IA,JA)

    integer  :: ifile, istep, t, tall
    integer  :: k, i, j, iq
    integer  :: ierr
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

    if ( BOUNDARY_UPDATE_DT <= 0.0_DP ) then
       write(*,*) 'xxx BOUNDARY_UPDATE_DT is necessary in real case preprocess'
       call PRC_MPIstop
    endif

    if ( FILETYPE_ORG == 'GrADS' ) then
       basename_mod = trim(BASENAME_ORG) ! namelist file name
    else
       if ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM ) then
          basename_mod = trim(BASENAME_ORG)//'_00000'
       else
          basename_mod = trim(BASENAME_ORG)
       endif
    endif

    call ParentAtmosSetup( FILETYPE_ORG,     & ![IN]
                           basename_mod,     & ![IN]
                           SERIAL_PROC_READ, & ![IN]
                           USE_FILE_DENSITY, & ![IN]
                           dims(:),          & ![OUT]
                           timelen           ) ![OUT]

    if ( timelen > 0 ) then
       NUMBER_OF_TSTEPS = timelen ! read from file
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Number of temporal data in each file : ', NUMBER_OF_TSTEPS

    do ifile = 1, NUMBER_OF_FILES

       if ( FILETYPE_ORG == 'GrADS' ) then
          if ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM ) then
             write(basename_mod,'(A,I5.5)') '_', ifile-1 ! only the number postfix
          else
             basename_mod = ''
          endif
       else
          if ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM ) then
             write(basename_mod,'(A,A,I5.5)') trim(BASENAME_ORG), '_', ifile-1
          else
             basename_mod = trim(BASENAME_ORG)
          endif
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** read external data from : ', trim(basename_mod)

       call ParentAtmosOpen( FILETYPE_ORG, & ![IN]
                             basename_mod, & ![IN]
                             dims(:)       ) ![IN]

       do istep = 1, NUMBER_OF_TSTEPS

          tall = NUMBER_OF_TSTEPS * (ifile-1) + istep ! consecutive time step (input)
          t    = tall - NUMBER_OF_SKIP_TSTEPS         ! time step (output)

          if ( t <= 0 ) then
             if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,I5,A,I6,A)') &
                        '*** [file,step,cons.] = [', ifile, ',', istep, ',', tall, '] ...skip.'
             cycle
          endif

          if ( t == 1 .OR. BASENAME_BOUNDARY /= '' ) then

             if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,I5,A,I6,A)') &
                        '*** [file,step,cons.] = [', ifile, ',', istep, ',', tall, ']'

             ! read prepared data
             call ParentAtmosInput( FILETYPE_ORG,     & ! [IN]
                                    basename_mod,     & ! [IN]
                                    dims(:),          & ! [IN]
                                    istep,            & ! [IN]
                                    SAME_MP_TYPE,     & ! [IN]
                                    DENS_in(:,:,:),   & ! [OUT]
                                    MOMZ_in(:,:,:),   & ! [OUT]
                                    MOMX_in(:,:,:),   & ! [OUT]
                                    MOMY_in(:,:,:),   & ! [OUT]
                                    RHOT_in(:,:,:),   & ! [OUT]
                                    QTRC_in(:,:,:,:), & ! [OUT]
                                    VELZ_in(:,:,:),   & ! [OUT]
                                    VELX_in(:,:,:),   & ! [OUT]
                                    VELY_in(:,:,:),   & ! [OUT]
                                    POTT_in(:,:,:)    ) ! [OUT]
          else
             if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,I5,A,I6,A)') &
                        '*** [file,step,cons.] = [', ifile, ',', istep, ',', tall, '] ...skip.'
          endif

          !--- store prognostic variables as initial
          if ( t == 1 ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** store initial state.'

             do j = 1, JA
             do i = 1, IA
             do k = 1, KA
                DENS(k,i,j) = DENS_in(k,i,j)
                MOMZ(k,i,j) = MOMZ_in(k,i,j)
                MOMX(k,i,j) = MOMX_in(k,i,j)
                MOMY(k,i,j) = MOMY_in(k,i,j)
                RHOT(k,i,j) = RHOT_in(k,i,j)
             enddo
             enddo
             enddo

             do iq = 1, QA
             do j  = 1, JA
             do i  = 1, IA
             do k  = 1, KA
                QTRC(k,i,j,iq) = QTRC_in(k,i,j,iq)
             enddo
             enddo
             enddo
             enddo

          endif

          !--- output boundary data
          if ( BASENAME_BOUNDARY /= '' ) then

             if ( t == 1 ) then
                if ( BOUNDARY_POSTFIX_TIMELABEL ) then
                   call TIME_gettimelabel( timelabel )
                   basename_out_mod = trim(BASENAME_BOUNDARY)//'_'//trim(timelabel)
                else
                   basename_out_mod = trim(BASENAME_BOUNDARY)
                endif

                call BoundaryAtmosSetup( basename_out_mod,   & ! [IN]
                                         BOUNDARY_TITLE,     & ! [IN]
                                         BOUNDARY_DTYPE,     & ! [IN]
                                         BOUNDARY_UPDATE_DT, & ! [IN]
                                         fid_atmos,          & ! [OUT]
                                         vid_atmos(:)        ) ! [OUT]
             endif

             call BoundaryAtmosOutput( DENS_in(:,:,:),     & ! [IN]
                                       VELZ_in(:,:,:),     & ! [IN]
                                       VELX_in(:,:,:),     & ! [IN]
                                       VELY_in(:,:,:),     & ! [IN]
                                       POTT_in(:,:,:),     & ! [IN]
                                       QTRC_in(:,:,:,:),   & ! [IN]
                                       fid_atmos,          & ! [IN]
                                       vid_atmos(:),       & ! [IN]
                                       BOUNDARY_UPDATE_DT, & ! [IN]
                                       t                   ) ! [IN]
          endif

       enddo ! istep loop
    enddo ! ifile loop

    return
  end subroutine REALINPUT_atmos

  !-----------------------------------------------------------------------------
  subroutine REALINPUT_surface
    use scale_const, only: &
         I_SW => CONST_I_SW, &
         I_LW => CONST_I_LW
    use scale_time, only: &
       TIME_gettimelabel
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
         OCEAN_SALT, &
         OCEAN_UVEL, &
         OCEAN_VVEL, &
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

    logical                  :: USE_FILE_LANDWATER   = .true.    ! use land water data from files
    real(RP)                 :: INIT_LANDWATER_RATIO = 0.5_RP    ! Ratio of land water to storage is constant, if USE_FILE_LANDWATER is ".false."
    real(RP)                 :: INIT_OCEAN_ALB_LW    = 0.04_RP   ! initial LW albedo on the ocean
    real(RP)                 :: INIT_OCEAN_ALB_SW    = 0.10_RP   ! initial SW albedo on the ocean
    real(RP)                 :: INIT_OCEAN_Z0W       = 1.0E-3_RP ! initial surface roughness on the ocean
    character(len=H_SHORT)   :: INTRP_LAND_TEMP      = 'off'
    character(len=H_SHORT)   :: INTRP_LAND_WATER     = 'off'
    character(len=H_SHORT)   :: INTRP_LAND_SFC_TEMP  = 'off'
    character(len=H_SHORT)   :: INTRP_OCEAN_TEMP     = 'off'
    character(len=H_SHORT)   :: INTRP_OCEAN_SFC_TEMP = 'off'
    integer                  :: INTRP_ITER_MAX       = 100
    character(len=H_SHORT)   :: SOILWATER_DS2VC      = 'limit'
    logical                  :: soilwater_DS2VC_flag           ! true: 'critical', false: 'limit'
    logical                  :: elevation_collection = .true.

    NAMELIST / PARAM_MKINIT_REAL_LAND / &
       NUMBER_OF_FILES,            &
       NUMBER_OF_TSTEPS,           &
       NUMBER_OF_SKIP_TSTEPS,      &
       FILETYPE_ORG,               &
       BASENAME_ORG,               &
       BASENAME_ADD_NUM,           &
       BASENAME_BOUNDARY,          &
       BOUNDARY_POSTFIX_TIMELABEL, &
       BOUNDARY_TITLE,             &
       BOUNDARY_UPDATE_DT,         &
       USE_FILE_LANDWATER,         &
       INIT_LANDWATER_RATIO,       &
       INTRP_LAND_TEMP,            &
       INTRP_LAND_WATER,           &
       INTRP_LAND_SFC_TEMP,        &
       INTRP_ITER_MAX,             &
       SOILWATER_DS2VC,            &
       ELEVATION_COLLECTION,       &
       SERIAL_PROC_READ

    NAMELIST / PARAM_MKINIT_REAL_OCEAN / &
       NUMBER_OF_FILES,            &
       NUMBER_OF_TSTEPS,           &
       NUMBER_OF_SKIP_TSTEPS,      &
       FILETYPE_ORG,               &
       BASENAME_ORG,               &
       BASENAME_ADD_NUM,           &
       BASENAME_BOUNDARY,          &
       BOUNDARY_POSTFIX_TIMELABEL, &
       BOUNDARY_TITLE,             &
       BOUNDARY_UPDATE_DT,         &
       INIT_OCEAN_ALB_LW,          &
       INIT_OCEAN_ALB_SW,          &
       INIT_OCEAN_Z0W,             &
       INTRP_OCEAN_TEMP,           &
       INTRP_OCEAN_SFC_TEMP,       &
       INTRP_ITER_MAX,             &
       SERIAL_PROC_READ

    character(len=H_LONG) :: FILETYPE_LAND
    character(len=H_LONG) :: FILETYPE_OCEAN
    character(len=H_LONG) :: BASENAME_LAND
    character(len=H_LONG) :: BASENAME_OCEAN
    character(len=5)      :: NUM               = ''

    ! land
    real(RP), allocatable :: LAND_TEMP_org      (:,:,:,:)
    real(RP), allocatable :: LAND_WATER_org     (:,:,:,:)
    real(RP), allocatable :: LAND_SFC_TEMP_org  (:,:,:)
    real(RP), allocatable :: LAND_SFC_albedo_org(:,:,:,:)

    ! urban
    real(RP) :: URBAN_TC_ORG(IA,JA)
    real(RP) :: URBAN_QC_ORG(IA,JA)
    real(RP) :: URBAN_UC_ORG(IA,JA)
    real(RP) :: URBAN_SFC_TEMP_ORG(IA,JA)
    real(RP) :: URBAN_SFC_albedo_ORG(IA,JA,2)

    ! ocean
    real(RP), allocatable :: OCEAN_TEMP_org      (:,:,:,:)
    real(RP), allocatable :: OCEAN_SFC_TEMP_org  (:,:,:)
    real(RP), allocatable :: OCEAN_SFC_albedo_org(:,:,:,:)
    real(RP), allocatable :: OCEAN_SFC_Z0_org    (:,:,:)

    integer :: NUMBER_OF_FILES_LAND        = 1
    integer :: NUMBER_OF_FILES_OCEAN       = 1
    integer :: NUMBER_OF_TSTEPS_LAND       = 1       ! num of time steps in one file
    integer :: NUMBER_OF_TSTEPS_OCEAN      = 1       ! num of time steps in one file
    integer :: NUMBER_OF_SKIP_TSTEPS_LAND  = 0       ! num of skipped first several data
    integer :: NUMBER_OF_SKIP_TSTEPS_OCEAN = 0       ! num of skipped first several data

    character(len=H_LONG) :: BASENAME_BOUNDARY_LAND           = ''
    character(len=H_LONG) :: BASENAME_BOUNDARY_OCEAN          = ''
    logical               :: BOUNDARY_POSTFIX_TIMELABEL_LAND  = .false.
    logical               :: BOUNDARY_POSTFIX_TIMELABEL_OCEAN = .false.
    character(len=H_LONG) :: BOUNDARY_TITLE_LAND              = 'SCALE-RM BOUNDARY CONDITION for REAL CASE'
    character(len=H_LONG) :: BOUNDARY_TITLE_OCEAN             = 'SCALE-RM BOUNDARY CONDITION for REAL CASE'
    real(DP)              :: BOUNDARY_UPDATE_DT_LAND          = 0.0_DP  ! inteval time of boudary data update [s]
    real(DP)              :: BOUNDARY_UPDATE_DT_OCEAN         = 0.0_DP  ! inteval time of boudary data update [s]

    integer :: mdlid_land, mdlid_ocean
    integer :: ldims(3), odims(2)

    integer :: totaltimesteps = 1
    integer :: timelen
    integer :: skip_steps
    integer :: ierr

    character(len=H_LONG) :: basename_out_mod
    character(len=19)     :: timelabel

    logical :: boundary_flag = .false.

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

    NUMBER_OF_FILES_LAND            = NUMBER_OF_FILES
    NUMBER_OF_TSTEPS_LAND           = NUMBER_OF_TSTEPS
    NUMBER_OF_SKIP_TSTEPS_LAND      = NUMBER_OF_SKIP_TSTEPS
    FILETYPE_LAND                   = FILETYPE_ORG
    BASENAME_BOUNDARY_LAND          = BASENAME_BOUNDARY
    BOUNDARY_POSTFIX_TIMELABEL_LAND = BOUNDARY_POSTFIX_TIMELABEL
    BOUNDARY_TITLE_LAND             = BOUNDARY_TITLE
    BOUNDARY_UPDATE_DT_LAND         = BOUNDARY_UPDATE_DT

    if ( FILETYPE_LAND .ne. "GrADS" .and. ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM ) ) then
       BASENAME_LAND = trim(BASENAME_ORG)//"_00000"
    else
       BASENAME_LAND = trim(BASENAME_ORG)
    endif

    select case( SOILWATER_DS2VC )
    case( 'critical' )
       SOILWATER_DS2VC_flag = .true.
    case('limit' )
       SOILWATER_DS2VC_flag = .false.
    case default
      write(*,*) 'xxx Unsupported SOILWATER_DS2CV TYPE:', trim(SOILWATER_DS2VC)
      call PRC_MPIstop
    end select

    serial_land = SERIAL_PROC_READ

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

    NUMBER_OF_FILES_OCEAN            = NUMBER_OF_FILES
    NUMBER_OF_TSTEPS_OCEAN           = NUMBER_OF_TSTEPS
    NUMBER_OF_SKIP_TSTEPS_OCEAN      = NUMBER_OF_SKIP_TSTEPS
    FILETYPE_OCEAN                   = FILETYPE_ORG
    BASENAME_BOUNDARY_OCEAN          = BASENAME_BOUNDARY
    BOUNDARY_POSTFIX_TIMELABEL_OCEAN = BOUNDARY_POSTFIX_TIMELABEL
    BOUNDARY_TITLE_OCEAN             = BOUNDARY_TITLE
    BOUNDARY_UPDATE_DT_OCEAN         = BOUNDARY_UPDATE_DT

    if ( FILETYPE_OCEAN .ne. "GrADS" .and. ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM ) ) then
       BASENAME_OCEAN = trim(BASENAME_ORG)//"_00000"
    else
       BASENAME_OCEAN = trim(BASENAME_ORG)
    endif

    serial_ocean = SERIAL_PROC_READ

    ! check land/ocean parameters
    if( NUMBER_OF_FILES_LAND            .NE.   NUMBER_OF_FILES_OCEAN            .OR. &
        NUMBER_OF_TSTEPS_LAND           .NE.   NUMBER_OF_TSTEPS_OCEAN           .OR. &
        NUMBER_OF_SKIP_TSTEPS_LAND      .NE.   NUMBER_OF_SKIP_TSTEPS_OCEAN      .OR. &
        BASENAME_BOUNDARY_LAND          .NE.   BASENAME_BOUNDARY_OCEAN          .OR. &
        BOUNDARY_POSTFIX_TIMELABEL_LAND .NEQV. BOUNDARY_POSTFIX_TIMELABEL_OCEAN .OR. &
        BOUNDARY_TITLE_LAND             .NE.   BOUNDARY_TITLE_OCEAN             .OR. &
        BOUNDARY_UPDATE_DT_LAND         .NE.   BOUNDARY_UPDATE_DT_OCEAN              ) then
       write(*,*) 'xxx Error: The following LAND/OCEAN parameters must be consistent due to technical problem:'
       write(*,*) '           NUMBER_OF_FILES, NUMBER_OF_TSTEPS, NUMBER_OF_SKIP_TSTEPS,'
       write(*,*) '           BASENAME_BOUNDARY, BOUNDARY_POSTFIX_TIMELABEL, BOUNDARY_TITLE, BOUNDARY_UPDATE_DT.'
       call PRC_MPIstop
    end if

    call ParentSurfaceSetup( ldims, odims,           & ![OUT]
                             mdlid_land,             & ![OUT]
                             mdlid_ocean,            & ![OUT]
                             timelen,                & ![OUT]
                             BASENAME_LAND,          & ![IN]
                             BASENAME_OCEAN,         & ![IN]
                             FILETYPE_LAND,          & ![IN]
                             FILETYPE_OCEAN,         & ![IN]
                             USE_FILE_LANDWATER,     & ![IN]
                             intrp_land_TEMP,        & ![IN]
                             intrp_land_water,       & ![IN]
                             intrp_land_sfc_TEMP,    & ![IN]
                             intrp_ocean_TEMP,       & ![IN]
                             intrp_ocean_sfc_TEMP    ) ![IN]

    if ( timelen > 0 ) then
       NUMBER_OF_TSTEPS = timelen ! read from file
    endif

    totaltimesteps = NUMBER_OF_FILES * NUMBER_OF_TSTEPS

    allocate( LAND_TEMP_ORG      (LKMAX,IA,JA,  1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )
    allocate( LAND_WATER_ORG     (LKMAX,IA,JA,  1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )
    allocate( LAND_SFC_TEMP_ORG  (      IA,JA,  1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )
    allocate( LAND_SFC_albedo_ORG(      IA,JA,2,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )

    allocate( OCEAN_TEMP_ORG      (OKMAX,IA,JA,  1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )
    allocate( OCEAN_SFC_TEMP_ORG  (      IA,JA,  1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )
    allocate( OCEAN_SFC_albedo_ORG(      IA,JA,2,1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )
    allocate( OCEAN_SFC_Z0_ORG    (      IA,JA,  1+NUMBER_OF_SKIP_TSTEPS:totaltimesteps) )

    if ( mdlid_ocean == iGrADS ) then
       BASENAME_ORG = ""
    endif

    if ( BASENAME_BOUNDARY /= '' ) then
       boundary_flag = .true.
    endif

    !--- read external file
    do n = 1, NUMBER_OF_FILES

       if ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM ) then
          write(NUM,'(I5.5)') n-1
          BASENAME_LAND  = trim(BASENAME_ORG)//"_"//NUM
          BASENAME_OCEAN = trim(BASENAME_ORG)//"_"//NUM
       else
          BASENAME_LAND  = trim(BASENAME_ORG)
          BASENAME_OCEAN = trim(BASENAME_ORG)
       endif

       if( IO_L ) write(IO_FID_LOG,*) ' '
       if( IO_L ) write(IO_FID_LOG,*) '+++ Target File Name (Land) : ', trim(BASENAME_LAND)
       if( IO_L ) write(IO_FID_LOG,*) '+++ Target File Name (Ocean): ', trim(BASENAME_OCEAN)
       if( IO_L ) write(IO_FID_LOG,*) '    Time Steps in One File  : ', NUMBER_OF_TSTEPS

       ns = NUMBER_OF_TSTEPS * (n - 1) + 1
       ne = ns + (NUMBER_OF_TSTEPS - 1)

       if ( ne <= NUMBER_OF_SKIP_TSTEPS ) then
          if( IO_L ) write(IO_FID_LOG,*) '    SKIP'
          cycle
       endif

       skip_steps = max(NUMBER_OF_SKIP_TSTEPS - ns + 1, 0)
       ns = max(ns, NUMBER_OF_SKIP_TSTEPS+1)

       ! read all prepared data
       call ParentSurfaceInput( LAND_TEMP_org       (:,:,:,ns:ne), &
                                LAND_WATER_org      (:,:,:,ns:ne), &
                                LAND_SFC_TEMP_org   (:,:,  ns:ne), &
                                LAND_SFC_albedo_org (:,:,:,ns:ne), &
                                URBAN_TC_org,         &
                                URBAN_QC_org,         &
                                URBAN_UC_org,         &
                                URBAN_SFC_TEMP_org,   &
                                URBAN_SFC_albedo_org, &
                                OCEAN_TEMP_org      (OKS,:,:,  ns:ne), &
                                OCEAN_SFC_TEMP_org  (    :,:,  ns:ne), &
                                OCEAN_SFC_albedo_org(    :,:,:,ns:ne), &
                                OCEAN_SFC_Z0_org    (    :,:,  ns:ne), &
                                BASENAME_LAND,           &
                                BASENAME_OCEAN,          &
                                mdlid_land, mdlid_ocean, &
                                ldims, odims,            &
                                USE_FILE_LANDWATER,      &
                                INIT_LANDWATER_RATIO,    &
                                INIT_OCEAN_ALB_LW,       &
                                INIT_OCEAN_ALB_SW,       &
                                INIT_OCEAN_Z0W,          &
                                INTRP_ITER_MAX,          &
                                SOILWATER_DS2VC_flag,    &
                                elevation_collection,    &
                                boundary_flag,           &
                                NUMBER_OF_TSTEPS,        &
                                skip_steps               )

       ! required one-step data only
       if( BASENAME_BOUNDARY == '' ) exit

    enddo


    !--- input initial data
    ns = NUMBER_OF_SKIP_TSTEPS + 1  ! skip first several data

    do j = 1, JA
    do i = 1, IA
       LAND_SFC_TEMP  (i,j)      = LAND_SFC_TEMP_org  (i,j,     ns)
       LAND_SFC_albedo(i,j,I_LW) = LAND_SFC_albedo_org(i,j,I_LW,ns)
       LAND_SFC_albedo(i,j,I_SW) = LAND_SFC_albedo_org(i,j,I_SW,ns)
       do k = 1, LKMAX
          LAND_TEMP (k,i,j) = LAND_TEMP_org (k,i,j,ns)
          LAND_WATER(k,i,j) = LAND_WATER_org(k,i,j,ns)
       enddo

       URBAN_SFC_TEMP  (i,j)      = URBAN_SFC_TEMP_org  (i,j)
       URBAN_SFC_albedo(i,j,I_LW) = URBAN_SFC_albedo_org(i,j,I_LW)
       URBAN_SFC_albedo(i,j,I_SW) = URBAN_SFC_albedo_org(i,j,I_SW)
       do k = UKS, UKE
          URBAN_TRL(k,i,j) = URBAN_SFC_TEMP_org(i,j)
          URBAN_TBL(k,i,j) = URBAN_SFC_TEMP_org(i,j)
          URBAN_TGL(k,i,j) = URBAN_SFC_TEMP_org(i,j)
       enddo
       URBAN_TC   (i,j) = URBAN_TC_org      (i,j)
       URBAN_QC   (i,j) = URBAN_QC_org      (i,j)
       URBAN_UC   (i,j) = URBAN_UC_org      (i,j)
       URBAN_TR   (i,j) = URBAN_SFC_TEMP_org(i,j)
       URBAN_TB   (i,j) = URBAN_SFC_TEMP_org(i,j)
       URBAN_TG   (i,j) = URBAN_SFC_TEMP_org(i,j)
       URBAN_RAINR(i,j) = 0.0_RP
       URBAN_RAINB(i,j) = 0.0_RP
       URBAN_RAING(i,j) = 0.0_RP
       URBAN_ROFF (i,j) = 0.0_RP

       do k = 1, OKMAX
          OCEAN_TEMP(k,i,j) = OCEAN_TEMP_ORG(OKS,i,j,ns)
          OCEAN_SALT(k,i,j) = 0.0_RP
          OCEAN_UVEL(k,i,j) = 0.0_RP
          OCEAN_VVEL(k,i,j) = 0.0_RP
       enddo
       OCEAN_SFC_TEMP  (i,j)      = OCEAN_SFC_TEMP_ORG  (i,j,     ns)
       OCEAN_SFC_albedo(i,j,I_LW) = OCEAN_SFC_albedo_ORG(i,j,I_LW,ns)
       OCEAN_SFC_albedo(i,j,I_SW) = OCEAN_SFC_albedo_ORG(i,j,I_SW,ns)
       OCEAN_SFC_Z0M   (i,j)      = OCEAN_SFC_Z0_ORG    (i,j,     ns)
       OCEAN_SFC_Z0H   (i,j)      = OCEAN_SFC_Z0_ORG    (i,j,     ns)
       OCEAN_SFC_Z0E   (i,j)      = OCEAN_SFC_Z0_ORG    (i,j,     ns)
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       ATMOS_PHY_SF_SFC_TEMP  (i,j)      = fact_ocean(i,j) * OCEAN_SFC_TEMP  (i,j) &
                                         + fact_land (i,j) * LAND_SFC_TEMP   (i,j) &
                                         + fact_urban(i,j) * URBAN_SFC_TEMP  (i,j)
       ATMOS_PHY_SF_SFC_albedo(i,j,I_LW) = fact_ocean(i,j) * OCEAN_SFC_albedo(i,j,I_LW) &
                                         + fact_land (i,j) * LAND_SFC_albedo (i,j,I_LW) &
                                         + fact_urban(i,j) * URBAN_SFC_albedo(i,j,I_LW)
       ATMOS_PHY_SF_SFC_albedo(i,j,I_SW) = fact_ocean(i,j) * OCEAN_SFC_albedo(i,j,I_SW) &
                                         + fact_land (i,j) * LAND_SFC_albedo (i,j,I_SW) &
                                         + fact_urban(i,j) * URBAN_SFC_albedo(i,j,I_SW)
       ATMOS_PHY_SF_SFC_Z0M   (i,j)      = OCEAN_SFC_Z0M(i,j)
       ATMOS_PHY_SF_SFC_Z0H   (i,j)      = OCEAN_SFC_Z0H(i,j)
       ATMOS_PHY_SF_SFC_Z0E   (i,j)      = OCEAN_SFC_Z0E(i,j)
    enddo
    enddo


    !--- output boundary data
    if( BASENAME_BOUNDARY /= '' ) then
       totaltimesteps = totaltimesteps - NUMBER_OF_SKIP_TSTEPS ! skip first several data
       if ( totaltimesteps > 1 ) then
          if ( BOUNDARY_UPDATE_DT <= 0.0_DP ) then
             write(*,*) 'xxx BOUNDARY_UPDATE_DT is necessary in real case preprocess'
             call PRC_MPIstop
          endif

          if ( BOUNDARY_POSTFIX_TIMELABEL ) then
             call TIME_gettimelabel( timelabel )
             basename_out_mod = trim(BASENAME_BOUNDARY)//'_'//trim(timelabel)
          else
             basename_out_mod = trim(BASENAME_BOUNDARY)
          endif

          call ParentSurfaceBoundary( LAND_TEMP_org       (:,:,:,  ns:ne), &
                                      LAND_WATER_org      (:,:,:,  ns:ne), &
                                      LAND_SFC_TEMP_org   (  :,:,  ns:ne), &
                                      LAND_SFC_albedo_org (  :,:,:,ns:ne), &
                                      OCEAN_TEMP_org      (:,:,:,  ns:ne), &
                                      OCEAN_SFC_TEMP_org  (  :,:,  ns:ne), &
                                      OCEAN_SFC_albedo_org(  :,:,:,ns:ne), &
                                      OCEAN_SFC_Z0_org    (  :,:,  ns:ne), &
                                      totaltimesteps,                    &
                                      BOUNDARY_UPDATE_DT,                &
                                      basename_out_mod,                  &
                                      BOUNDARY_TITLE                     )

       endif
    endif

    deallocate( LAND_TEMP_org        )
    deallocate( LAND_WATER_org       )
    deallocate( LAND_SFC_TEMP_org    )
    deallocate( LAND_SFC_albedo_org  )
    deallocate( OCEAN_TEMP_org       )
    deallocate( OCEAN_SFC_TEMP_org   )
    deallocate( OCEAN_SFC_albedo_org )
    deallocate( OCEAN_SFC_Z0_org     )

    return
  end subroutine REALINPUT_surface


  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtmosSetup( &
       inputtype,           &
       basename,            &
       serial_in,           &
       use_file_density_in, &
       dims,                &
       timelen              )
    use mod_realinput_scale, only: &
       ParentAtmosSetupSCALE
    use mod_realinput_wrfarw, only: &
       ParentAtmosSetupWRFARW
!!$    use mod_realinput_nicam, only: &
!!$       ParentAtmosSetupNICAM
    use mod_realinput_grads, only: &
       ParentAtmosSetupGrADS
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    character(len=*), intent(in)  :: inputtype
    character(len=*), intent(in)  :: basename
    logical,          intent(in)  :: serial_in           ! read by a serial process
    logical,          intent(in)  :: use_file_density_in ! use density data from files
    integer,          intent(out) :: dims(6)
    integer,          intent(out) :: timelen
    !---------------------------------------------------------------------------

    serial_atmos = serial_in
    if ( serial_atmos ) then
       if( PRC_IsMaster ) then
          read_by_myproc_atmos = .true.
       else
          read_by_myproc_atmos = .false.
       endif
    else
       read_by_myproc_atmos = .true.
    endif

    select case(inputtype)
    case('SCALE-RM')

       serial_atmos         = .false. ! force false
       read_by_myproc_atmos = .true.

       call ParentAtmosSetupSCALE( dims(:) ) ! [OUT]
       timelen = -1

       use_file_density = use_file_density_in
       temp2pott        = .false.
       update_coord     = .false.
       apply_rotate_uv  = .false.

    case('GrADS')

       if ( read_by_myproc_atmos ) then
          call ParentAtmosSetupGrADS ( dims(:), & ! [OUT]
                                       basename ) ! [IN]
       endif
       timelen = -1

       use_file_density = use_file_density_in
       temp2pott        = .true.
       update_coord     = .true.
       apply_rotate_uv  = .true.

    case('WRF-ARW')

       if ( read_by_myproc_atmos ) then
          call ParentAtmosSetupWRFARW( dims(:), & ! [OUT]
                                       timelen, & ! [OUT]
                                       basename ) ! [IN]
       endif

       use_file_density = .false.
       temp2pott        = .true.
       update_coord     = .true.
       apply_rotate_uv  = .true.

!!$    case('NICAM-NETCDF')
!!$
!!$       if ( read_by_myproc_atmos ) then
!!$          call ParentAtmosSetupNICAM ( dims(:), & ! [OUT]
!!$                                       timelen, & ! [OUT]
!!$                                       basename ) ! [IN]
!!$       endif
!!$
!!$       use_file_density = .false.
!!$       temp2pott        = .true.
!!$       update_coord     = .false.
!!$       apply_rotate_uv  = .true.
!!$
    case default

       write(*,*) 'xxx Unsupported type of input data : ', trim(inputtype)
       call PRC_MPIstop

    end select

    if ( serial_atmos ) then
       call COMM_bcast( dims(:), 6 )
       call COMM_bcast( timelen )
    endif

    allocate( LON_org (            dims(2), dims(3)     ) )
    allocate( LAT_org (            dims(2), dims(3)     ) )
    allocate( CZ_org  ( dims(1)+2, dims(2), dims(3)     ) )

    allocate( W_org   ( dims(1)+2, dims(2), dims(3)     ) )
    allocate( U_org   ( dims(1)+2, dims(2), dims(3)     ) )
    allocate( V_org   ( dims(1)+2, dims(2), dims(3)     ) )
    allocate( POTT_org( dims(1)+2, dims(2), dims(3)     ) )
    allocate( TEMP_org( dims(1)+2, dims(2), dims(3)     ) )
    allocate( PRES_org( dims(1)+2, dims(2), dims(3)     ) )
    allocate( DENS_org( dims(1)+2, dims(2), dims(3)     ) )
    allocate( QTRC_org( dims(1)+2, dims(2), dims(3), QA ) )

    allocate( QV_org  ( dims(1)+2, dims(2), dims(3)        ) )
    allocate( QHYD_org( dims(1)+2, dims(2), dims(3), N_HYD ) )
    allocate( QNUM_org( dims(1)+2, dims(2), dims(3), N_HYD ) )

    if( IO_L ) write(IO_FID_LOG,*) '*** Horizontal Interpolation Level: ', COMM_CARTESC_NEST_INTERP_LEVEL
    itp_nh = COMM_CARTESC_NEST_INTERP_LEVEL
    itp_nv = 2

    allocate( igrd (          IA,JA,itp_nh) )
    allocate( jgrd (          IA,JA,itp_nh) )
    allocate( hfact(          IA,JA,itp_nh) )
    allocate( kgrd (KA,itp_nv,IA,JA,itp_nh) )
    allocate( vfact(KA,itp_nv,IA,JA,itp_nh) )

    return
  end subroutine ParentAtmosSetup

  !-----------------------------------------------------------------------------
  !> Atmosphere Data Open
  subroutine ParentAtmosOpen( &
       inputtype, &
       basename,  &
       dims       )
    use mod_realinput_scale, only: &
       ParentAtmosOpenSCALE
    use mod_realinput_wrfarw, only: &
       ParentAtmosOpenWRFARW
!!$    use mod_realinput_nicam, only: &
!!$       ParentAtmosOpenNICAM
    use mod_realinput_grads, only: &
       ParentAtmosOpenGrADS
    implicit none

    character(len=*), intent(in)  :: inputtype
    character(len=*), intent(in)  :: basename
    integer,          intent(in)  :: dims(6)
    !---------------------------------------------------------------------------

    if ( read_by_myproc_atmos ) then

       select case(inputtype)
       case('SCALE-RM')
          call ParentAtmosOpenSCALE( LON_org(:,:),   & ! [OUT]
                                     LAT_org(:,:),   & ! [OUT]
                                     CZ_org (:,:,:), & ! [OUT]
                                     basename,       & ! [IN]
                                     dims   (:)      ) ! [IN]
       case('GrADS')
          call ParentAtmosOpenGrADS
       case('WRF-ARW')
          call ParentAtmosOpenWRFARW
       case('NETCDF')
          call ParentAtmosOpenSCALE( LON_org(:,:),   & ! [OUT]
                                     LAT_org(:,:),   & ! [OUT]
                                     CZ_org (:,:,:), & ! [OUT]
                                     basename,       & ! [IN]
                                     dims   (:)      ) ! [IN]
       end select

    endif

    return
  end subroutine ParentAtmosOpen

  !-----------------------------------------------------------------------------
  !> Atmosphere Data Read
  subroutine ParentAtmosInput( &
       inputtype,     &
       basename,      &
       dims,          &
       istep,         &
       same_mptype,   &
       DENS,          &
       MOMZ,          &
       MOMX,          &
       MOMY,          &
       RHOT,          &
       QTRC,          &
       VELZ,          &
       VELX,          &
       VELY,          &
       POTT           )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_grid_cartesC_metric, only: &
       ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       I_QV, &
       QLS, &
       QLE
    use scale_atmos_thermodyn, only: &
       THERMODYN_qdry           => ATMOS_THERMODYN_qdry, &
       THERMODYN_r              => ATMOS_THERMODYN_r,    &
       THERMODYN_cp             => ATMOS_THERMODYN_cp,   &
       THERMODYN_temp_pres2pott => ATMOS_THERMODYN_temp_pres2pott
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_interp, only: &
       INTRP_domain_compatibility, &
       INTRP_factor3d,             &
       INTRP_interp3d
    use mod_realinput_scale, only: &
       ParentAtmosInputSCALE
    use mod_realinput_wrfarw, only: &
       ParentAtmosInputWRFARW
!!$    use mod_realinput_nicam, only: &
!!$       ParentAtmosInputNICAM
    use mod_realinput_grads, only: &
       ParentAtmosInputGrADS
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    use scale_atmos_grid_cartesC_real, only: &
       CZ => ATMOS_GRID_CARTESC_REAL_CZ
    implicit none

    character(len=*), intent(in)  :: inputtype
    character(len=*), intent(in)  :: basename
    integer,          intent(in)  :: dims(6)
    integer,          intent(in)  :: istep
    logical,          intent(in)  :: same_mptype   ! Is microphysics type same between outer and inner model
    real(RP),         intent(out) :: DENS(KA,IA,JA)
    real(RP),         intent(out) :: MOMZ(KA,IA,JA)
    real(RP),         intent(out) :: MOMX(KA,IA,JA)
    real(RP),         intent(out) :: MOMY(KA,IA,JA)
    real(RP),         intent(out) :: RHOT(KA,IA,JA)
    real(RP),         intent(out) :: QTRC(KA,IA,JA,QA)
    real(RP),         intent(out) :: VELZ(KA,IA,JA)
    real(RP),         intent(out) :: VELX(KA,IA,JA)
    real(RP),         intent(out) :: VELY(KA,IA,JA)
    real(RP),         intent(out) :: POTT(KA,IA,JA)

    real(RP) :: PRES (KA,IA,JA)
    real(RP) :: TEMP (KA,IA,JA)
    real(RP) :: W    (KA,IA,JA)
    real(RP) :: U    (KA,IA,JA)
    real(RP) :: V    (KA,IA,JA)
    real(RP) :: QV   (KA,IA,JA)
    real(RP) :: QC   (KA,IA,JA)
    real(RP) :: u_on_map, v_on_map

    real(RP) :: qdry, Rtot, CPtot

    logical, save :: first = .true.

    logical :: same_mptype_ = .false.

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( read_by_myproc_atmos ) then
       call PROF_rapstart('___AtmosInput',3)

       select case(inputtype)
       case('SCALE-RM')
          call ParentAtmosInputSCALE ( W_org   (:,:,:),   & ! [OUT]
                                       U_org   (:,:,:),   & ! [OUT]
                                       V_org   (:,:,:),   & ! [OUT]
                                       PRES_org(:,:,:),   & ! [OUT]
                                       DENS_org(:,:,:),   & ! [OUT]
                                       POTT_org(:,:,:),   & ! [OUT]
                                       QV_org  (:,:,:),   & ! [OUT]
                                       QHYD_org(:,:,:,:), & ! [OUT]
                                       QNUM_org(:,:,:,:), & ! [OUT]
                                       QTRC_org(:,:,:,:), & ! [OUT]
                                       CZ_org  (:,:,:),   & ! [IN]
                                       basename,          & ! [IN]
                                       same_mptype,       & ! [IN]
                                       dims(:),           & ! [IN]
                                       istep              ) ! [IN]
          same_mptype_ = same_mptype
       case('GrADS')
          call ParentAtmosInputGrADS ( W_org   (:,:,:),   & ! [OUT]
                                       U_org   (:,:,:),   & ! [OUT]
                                       V_org   (:,:,:),   & ! [OUT]
                                       PRES_org(:,:,:),   & ! [OUT]
                                       DENS_org(:,:,:),   & ! [OUT]
                                       TEMP_org(:,:,:),   & ! [OUT]
                                       QV_org  (:,:,:),   & ! [OUT]
                                       QHYD_org(:,:,:,:), & ! [OUT]
                                       LON_org (:,:),     & ! [OUT]
                                       LAT_org (:,:),     & ! [OUT]
                                       CZ_org  (:,:,:),   & ! [OUT]
                                       basename,          & ! [IN]
                                       dims(:),           & ! [IN]
                                       istep              ) ! [IN]
          same_mptype_ = .false.
          QNUM_org(:,:,:,:) = 0.0_RP
       case('WRF-ARW')
          call ParentAtmosInputWRFARW( W_org   (:,:,:),   & ! [OUT]
                                       U_org   (:,:,:),   & ! [OUT]
                                       V_org   (:,:,:),   & ! [OUT]
                                       PRES_org(:,:,:),   & ! [OUT]
                                       TEMP_org(:,:,:),   & ! [OUT]
                                       QV_org  (:,:,:),   & ! [OUT]
                                       QHYD_org(:,:,:,:), & ! [OUT]
                                       QNUM_org(:,:,:,:), & ! [OUT]
                                       LON_org (:,:),     & ! [OUT]
                                       LAT_org (:,:),     & ! [OUT]
                                       CZ_org  (:,:,:),   & ! [OUT]
                                       basename,          & ! [IN]
                                       dims(:),           & ! [IN]
                                       istep              ) ! [IN]
          same_mptype_ = .false.
          DENS_org(:,:,:) = 0.0_RP
!!$       case('NETCDF')
!!$          call ParentAtmosInputNICAM ( W_org   (:,:,:),   & ! [OUT]
!!$                                       U_org   (:,:,:),   & ! [OUT]
!!$                                       V_org   (:,:,:),   & ! [OUT]
!!$                                       PRES_org(:,:,:),   & ! [OUT]
!!$                                       TEMP_org(:,:,:),   & ! [OUT]
!!$                                       QTRC_org(:,:,:,:), & ! [OUT]
!!$                                       basename,          & ! [IN]
!!$                                       dims(:),           & ! [IN]
!!$                                       istep              ) ! [IN]
!!$          DENS_org(:,:,:) = 0.0_RP
       end select

       if ( .not. same_mptype_ ) then
          call ATMOS_PHY_MP_driver_qhyd2qtrc( dims(1)+2, 1, dims(1)+2, dims(2), 1, dims(2), dims(3), 1, dims(3), &
                                              QV_org(:,:,:), QHYD_org(:,:,:,:), & ! [IN]
                                              QTRC_org(:,:,:,QS_MP:QE_MP),      & ! [OUT]
                                              QNUM=QNUM_org(:,:,:,:)            ) ! [IN]
       end if

       if ( temp2pott ) then
          do j = 1, dims(3)
          do i = 1, dims(2)
          do k = 1, dims(1)+2
             call THERMODYN_qdry( QA, QTRC_org(k,i,j,:), TRACER_MASS(:), qdry )
             call THERMODYN_r   ( QA, QTRC_org(k,i,j,:), TRACER_R(:), qdry, Rtot )
             call THERMODYN_cp  ( QA, QTRC_org(k,i,j,:), TRACER_CP(:), qdry, CPtot )
             call THERMODYN_temp_pres2pott( TEMP_org(k,i,j), PRES_org(k,i,j), CPtot, Rtot, & ! [IN]
                                            POTT_org(k,i,j)                                ) ! [OUT]
          enddo
          enddo
          enddo
       endif

       call PROF_rapend  ('___AtmosInput',3)
    endif ! read by this process?

    if ( serial_atmos ) then
       call PROF_rapstart('___AtmosBcast',3)

       if ( first .OR. update_coord ) then
          call COMM_bcast( LON_org,            dims(2), dims(3) )
          call COMM_bcast( LAT_org,            dims(2), dims(3) )
          call COMM_bcast( CZ_org,  dims(1)+2, dims(2), dims(3) )
       endif

       call COMM_bcast( W_org   , dims(1)+2, dims(2), dims(3) )
       call COMM_bcast( U_org   , dims(1)+2, dims(2), dims(3) )
       call COMM_bcast( V_org   , dims(1)+2, dims(2), dims(3) )
       call COMM_bcast( POTT_org, dims(1)+2, dims(2), dims(3) )
       call COMM_bcast( PRES_org, dims(1)+2, dims(2), dims(3) )
       call COMM_bcast( DENS_org, dims(1)+2, dims(2), dims(3) )
       call COMM_bcast( QTRC_org, dims(1)+2, dims(2), dims(3), QA )

       call PROF_rapend  ('___AtmosBcast',3)
    endif

    do iq = 1, QA
    do j  = 1, dims(3)
    do i  = 1, dims(2)
    do k  = 1, dims(1)+2
       QTRC_org(k,i,j,iq) = max( QTRC_org(k,i,j,iq), 0.0_RP )
    enddo
    enddo
    enddo
    enddo

    ! interpolation
    call PROF_rapstart('___AtmosInterp',3)

    if ( first .OR. update_coord ) then
       first = .false.

       call INTRP_domain_compatibility( LON_org(:,:),      & ! [IN]
                                        LAT_org(:,:),      & ! [IN]
                                        CZ_org (:,:,:),    & ! [IN]
                                        LON    (:,:),      & ! [IN]
                                        LAT    (:,:),      & ! [IN]
                                        CZ     (KS:KE,:,:) ) ! [IN]

       ! full level
       call INTRP_factor3d( itp_nh,                  & ! [IN]
                            dims(1)+2, 1, dims(1)+2, & ! [IN]
                            dims(2), dims(3),        & ! [IN]
                            LON_org(:,:),            & ! [IN]
                            LAT_org(:,:),            & ! [IN]
                            CZ_org (:,:,:),          & ! [IN]
                            KA, KS, KE,              & ! [IN]
                            IA, JA,                  & ! [IN]
                            LON    (:,:),            & ! [IN]
                            LAT    (:,:),            & ! [IN]
                            CZ     (:,:,:),          & ! [IN]
                            igrd   (    :,:,:),      & ! [OUT]
                            jgrd   (    :,:,:),      & ! [OUT]
                            hfact  (    :,:,:),      & ! [OUT]
                            kgrd   (:,:,:,:,:),      & ! [OUT]
                            vfact  (:,:,:,:,:)       ) ! [OUT]
    endif

    call INTRP_interp3d( itp_nh,                      & ! [IN]
                         dims(1)+2, dims(2), dims(3), & ! [IN]
                         KA, KS, KE,                  & ! [IN]
                         IA, JA,                      & ! [IN]
                         igrd    (    :,:,:),         & ! [IN]
                         jgrd    (    :,:,:),         & ! [IN]
                         hfact   (    :,:,:),         & ! [IN]
                         kgrd    (:,:,:,:,:),         & ! [IN]
                         vfact   (:,:,:,:,:),         & ! [IN]
                         W_org   (:,:,:),             & ! [IN]
                         W       (:,:,:)              ) ! [OUT]

    call INTRP_interp3d( itp_nh,                      & ! [IN]
                         dims(1)+2, dims(2), dims(3), & ! [IN]
                         KA, KS, KE,                  & ! [IN]
                         IA, JA,                      & ! [IN]
                         igrd    (    :,:,:),         & ! [IN]
                         jgrd    (    :,:,:),         & ! [IN]
                         hfact   (    :,:,:),         & ! [IN]
                         kgrd    (:,:,:,:,:),         & ! [IN]
                         vfact   (:,:,:,:,:),         & ! [IN]
                         U_org   (:,:,:),             & ! [IN]
                         U       (:,:,:)              ) ! [OUT]

    call INTRP_interp3d( itp_nh,                      & ! [IN]
                         dims(1)+2, dims(2), dims(3), & ! [IN]
                         KA, KS, KE,                  & ! [IN]
                         IA, JA,                      & ! [IN]
                         igrd    (    :,:,:),         & ! [IN]
                         jgrd    (    :,:,:),         & ! [IN]
                         hfact   (    :,:,:),         & ! [IN]
                         kgrd    (:,:,:,:,:),         & ! [IN]
                         vfact   (:,:,:,:,:),         & ! [IN]
                         V_org   (:,:,:),             & ! [IN]
                         V       (:,:,:)              ) ! [OUT]

    if ( apply_rotate_uv ) then ! rotation from latlon field to map-projected field
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          u_on_map =  U(k,i,j) * ROTC(i,j,1) + V(k,i,j) * ROTC(i,j,2)
          v_on_map = -U(k,i,j) * ROTC(i,j,2) + V(k,i,j) * ROTC(i,j,1)

          U(k,i,j) = u_on_map
          V(k,i,j) = v_on_map
       enddo
       enddo
       enddo
    endif

    ! from scalar point to staggered point
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE-1
       VELZ(k,i,j) = 0.5_RP * ( W(k+1,i,j) + W(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA-1
    do k = KS, KE
       VELX(k,i,j) = 0.5_RP * ( U(k,i+1,j) + U(k,i,j) )
    enddo
    enddo
    enddo

    i = IA
    do j = 1, JA
    do k = KS, KE
       VELX(k,i,j) = U(k,i,j)
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = KS, KE
       VELY(k,i,j) = 0.5_RP * ( V(k,i,j+1) + V(k,i,j) )
    enddo
    enddo
    enddo

    j = JA
    do i = 1, IA
    do k = KS, KE
       VELY(k,i,j) = V(k,i,j)
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       VELZ(   1:KS-1,i,j) = 0.0_RP
       VELZ(KE  :KA  ,i,j) = 0.0_RP
       VELX(   1:KS-1,i,j) = 0.0_RP
       VELX(KE+1:KA  ,i,j) = 0.0_RP
       VELY(   1:KS-1,i,j) = 0.0_RP
       VELY(KE+1:KA  ,i,j) = 0.0_RP
    enddo
    enddo

    call COMM_vars8( VELZ(:,:,:), 1 )
    call COMM_vars8( VELX(:,:,:), 2 )
    call COMM_vars8( VELY(:,:,:), 3 )
    call COMM_wait ( VELZ(:,:,:), 1, .false. )
    call COMM_wait ( VELX(:,:,:), 2, .false. )
    call COMM_wait ( VELY(:,:,:), 3, .false. )

    call INTRP_interp3d( itp_nh,                      & ! [IN]
                         dims(1)+2, dims(2), dims(3), & ! [IN]
                         KA, KS, KE,                  & ! [IN]
                         IA, JA,                      & ! [IN]
                         igrd    (    :,:,:),         & ! [IN]
                         jgrd    (    :,:,:),         & ! [IN]
                         hfact   (    :,:,:),         & ! [IN]
                         kgrd    (:,:,:,:,:),         & ! [IN]
                         vfact   (:,:,:,:,:),         & ! [IN]
                         POTT_org(:,:,:),             & ! [IN]
                         POTT    (:,:,:)              ) ! [OUT]

    do j = 1, JA
    do i = 1, IA
       POTT(   1:KS-1,i,j) = 0.0_RP
       POTT(KE+1:KA  ,i,j) = 0.0_RP
    enddo
    enddo

    do iq = 1, QA
       call INTRP_interp3d( itp_nh,                      & ! [IN]
                            dims(1)+2, dims(2), dims(3), & ! [IN]
                            KA, KS, KE,                  & ! [IN]
                            IA, JA,                      & ! [IN]
                            igrd    (    :,:,:),         & ! [IN]
                            jgrd    (    :,:,:),         & ! [IN]
                            hfact   (    :,:,:),         & ! [IN]
                            kgrd    (:,:,:,:,:),         & ! [IN]
                            vfact   (:,:,:,:,:),         & ! [IN]
                            QTRC_org(:,:,:,iq),          & ! [IN]
                            QTRC    (:,:,:,iq)           ) ! [OUT]

       do j = 1, JA
       do i = 1, IA
          QTRC(   1:KS-1,i,j,iq) = 0.0_RP
          QTRC(KE+1:KA  ,i,j,iq) = 0.0_RP
       enddo
       enddo
    enddo

    if ( use_file_density ) then
       call INTRP_interp3d( itp_nh,                      & ! [IN]
                            dims(1)+2, dims(2), dims(3), & ! [IN]
                            KA, KS, KE,                  & ! [IN]
                            IA, JA,                      & ! [IN]
                            igrd    (    :,:,:),         & ! [IN]
                            jgrd    (    :,:,:),         & ! [IN]
                            hfact   (    :,:,:),         & ! [IN]
                            kgrd    (:,:,:,:,:),         & ! [IN]
                            vfact   (:,:,:,:,:),         & ! [IN]
                            DENS_org(:,:,:),             & ! [IN]
                            DENS    (:,:,:),             & ! [OUT]
                            logwgt = .true.              ) ! [IN]
    else
       call INTRP_interp3d( itp_nh,                      & ! [IN]
                            dims(1)+2, dims(2), dims(3), & ! [IN]
                            KA, KS, KE,                  & ! [IN]
                            IA, JA,                      & ! [IN]
                            igrd    (    :,:,:),         & ! [IN]
                            jgrd    (    :,:,:),         & ! [IN]
                            hfact   (    :,:,:),         & ! [IN]
                            kgrd    (:,:,:,:,:),         & ! [IN]
                            vfact   (:,:,:,:,:),         & ! [IN]
                            PRES_org(:,:,:),             & ! [IN]
                            PRES    (:,:,:),             & ! [OUT]
                            logwgt = .true.              ) ! [IN]

       QC(:,:,:) = 0.0_RP
       if ( ATMOS_HYDROMETEOR_dry ) then
          QV(:,:,:) = 0.0_RP
       else
          QV(:,:,:) = QTRC(:,:,:,I_QV)
          do iq = QLS, QLE
             QC(:,:,:) = QC(:,:,:) + QTRC(:,:,:,iq)
          enddo
       end if

       ! make density & pressure profile in moist condition
       call HYDROSTATIC_buildrho_real( KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                                       POTT(:,:,:), QV(:,:,:), QC(:,:,:), & ! [IN]
                                       CZ(:,:,:),                         & ! [IN]
                                       PRES(:,:,:),                       & ! [INOUT]
                                       DENS(:,:,:), TEMP(:,:,:)           ) ! [OUT]

       call COMM_vars8( DENS(:,:,:), 1 )
       call COMM_wait ( DENS(:,:,:), 1 )
    endif

    do j = 1, JA
    do i = 1, IA
       DENS(   1:KS-1,i,j) = 0.0_RP
       DENS(KE+1:KA  ,i,j) = 0.0_RP
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE-1
       MOMZ(k,i,j) = VELZ(k,i,j) * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA-1
    do k = KS, KE
       MOMX(k,i,j) = VELX(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    i = IA
    do j = 1, JA
    do k = KS, KE
       MOMX(k,i,j) = VELX(k,i,j) * DENS(k,i,j)
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = KS, KE
       MOMY(k,i,j) = VELY(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    j = JA
    do i = 1, IA
    do k = KS, KE
       MOMY(k,i,j) = VELY(k,i,j) * DENS(k,i,j)
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RHOT(k,i,j) = POTT(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       MOMZ(   1:KS-1,i,j) = 0.0_RP
       MOMZ(KE  :KA  ,i,j) = 0.0_RP
       MOMX(   1:KS-1,i,j) = 0.0_RP
       MOMX(KE+1:KA  ,i,j) = 0.0_RP
       MOMY(   1:KS-1,i,j) = 0.0_RP
       MOMY(KE+1:KA  ,i,j) = 0.0_RP
    enddo
    enddo

    call COMM_vars8( MOMZ(:,:,:), 1 )
    call COMM_vars8( MOMX(:,:,:), 2 )
    call COMM_vars8( MOMY(:,:,:), 3 )
    call COMM_wait ( MOMZ(:,:,:), 1, .false. )
    call COMM_wait ( MOMX(:,:,:), 2, .false. )
    call COMM_wait ( MOMY(:,:,:), 3, .false. )

    call PROF_rapend  ('___AtmosInterp',3)

    return
  end subroutine ParentAtmosInput

  !-----------------------------------------------------------------------------
  !> Boundary Data Write
  subroutine BoundaryAtmosSetup( &
       basename,  &
       title,     &
       datatype,  &
       timeintv,  &
       fid,       &
       vid        )
    use scale_file_cartesC, only: &
       FILE_CARTESC_create, &
       FILE_CARTESC_def_var, &
       FILE_CARTESC_enddef
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none

    character(len=*), intent(in)  :: basename
    character(len=*), intent(in)  :: title
    character(len=*), intent(in)  :: datatype
    real(DP),         intent(in)  :: timeintv
    integer,          intent(out) :: fid
    integer,          intent(out) :: vid(5+QA)

    integer :: iq
    !---------------------------------------------------------------------------

    call FILE_CARTESC_create( basename, title, datatype, fid, date=NOWDATE )

    call FILE_CARTESC_def_var( fid, &
         'DENS', 'Reference Density', 'kg/m3', 'ZXYT',  datatype, & ! [IN]
         vid(1),                                                  & ! [OUT]
         timeintv=timeintv                                        ) ! [IN]
    call FILE_CARTESC_def_var( fid, &
         'VELZ', 'Reference VELZ',    'm/s',   'ZHXYT', datatype, & ! [IN]
         vid(2),                                                  & ! [OUT]
         timeintv=timeintv                                        ) ! [IN]
    call FILE_CARTESC_def_var( fid, &
         'VELX', 'Reference VELX',    'm/s',   'ZXHYT', datatype, & ! [IN]
         vid(3),                                                  & ! [OUT]
         timeintv=timeintv                                        ) ! [IN]
    call FILE_CARTESC_def_var( fid, &
         'VELY', 'Reference VELY',    'm/s',   'ZXYHT', datatype, & ! [IN]
         vid(4),                                                  & ! [OUT]
         timeintv=timeintv                                        ) ! [IN]
    call FILE_CARTESC_def_var( fid, &
         'POTT', 'Reference PT',      'K',     'ZXYT',  datatype, & ! [IN]
         vid(5),                                                  & ! [OUT]
         timeintv=timeintv                                        ) ! [IN]

    do iq = QS_MP, QE_MP
       call FILE_CARTESC_def_var( fid,                               & ! [IN]
            TRACER_NAME(iq), 'Reference '//TRACER_NAME(iq), 'kg/kg', & ! [IN]
            'ZXYT', datatype,                                        & ! [IN]
            vid(5+iq),                                               & ! [OUT]
            timeintv = timeintv                                      ) ! [IN]
    enddo

    call FILE_CARTESC_enddef( fid )

    return
  end subroutine BoundaryAtmosSetup

  !-----------------------------------------------------------------------------
  !> Boundary Data Write
  subroutine BoundaryAtmosOutput( &
       DENS,     &
       VELZ,     &
       VELX,     &
       VELY,     &
       POTT,     &
       QTRC,     &
       fid,      &
       vid,      &
       timeintv, &
       istep     )
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: VELZ(KA,IA,JA)
    real(RP), intent(in)  :: VELX(KA,IA,JA)
    real(RP), intent(in)  :: VELY(KA,IA,JA)
    real(RP), intent(in)  :: POTT(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
    integer,  intent(in)  :: fid
    integer,  intent(in)  :: vid(5+QA)
    real(DP), intent(in)  :: timeintv
    integer,  intent(in)  :: istep

    real(RP) :: work(KA,IA,JA,1)

    real(DP) :: timeofs
    integer  :: iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('___AtmosOutput',3)

    timeofs = real(istep-1,kind=DP) * timeintv

!OCL XFILL
    work(:,:,:,1) = DENS(:,:,:)
    call FILE_CARTESC_write_var( fid, vid(1), work(:,:,:,:), 'DENS', 'ZXYT', timeintv, timeofs=timeofs )
!OCL XFILL
    work(:,:,:,1) = VELZ(:,:,:)
    call FILE_CARTESC_write_var( fid, vid(2), work(:,:,:,:), 'VELZ', 'ZHXYT', timeintv, timeofs=timeofs )
!OCL XFILL
    work(:,:,:,1) = VELX(:,:,:)
    call FILE_CARTESC_write_var( fid, vid(3), work(:,:,:,:), 'VELX', 'ZXHYT', timeintv, timeofs=timeofs )
!OCL XFILL
    work(:,:,:,1) = VELY(:,:,:)
    call FILE_CARTESC_write_var( fid, vid(4), work(:,:,:,:), 'VELY', 'ZXYHT', timeintv, timeofs=timeofs )
!OCL XFILL
    work(:,:,:,1) = POTT(:,:,:)
    call FILE_CARTESC_write_var( fid, vid(5), work(:,:,:,:), 'POTT', 'ZXYT', timeintv, timeofs=timeofs )

    do iq = QS_MP, QE_MP
       call FILE_CARTESC_write_var( fid, vid(5+iq),QTRC(:,:,:,iq:iq), TRACER_NAME(iq), &
                              'ZXYT', timeintv, timeofs=timeofs                  )
    enddo

    call PROF_rapend  ('___AtmosOutput',3)

    return
  end subroutine BoundaryAtmosOutput

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
    use mod_realinput_scale, only: &
         ParentLandSetupSCALE, &
         ParentOceanSetupSCALE
    use mod_realinput_wrfarw, only: &
         ParentLandSetupWRFARW, &
         ParentOceanSetupWRFARW
!!$    use mod_realinput_nicam, only: &
!!$         ParentLandSetupNICAM, &
!!$         ParentOceanSetupNICAM
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
       use_waterratio = .false.

!!$    case('NICAM-NETCDF')
!!$
!!$       lmdlid = iNICAM
!!$       if ( do_read_land ) call ParentLandSetupNICAM( ldims,        & ! (out)
!!$                                                      basename_land ) ! (in)
!!$       use_waterratio = .false.
!!$
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
!    case( iSCALE, iWRFARW, iNICAM )
    case( iSCALE, iWRFARW )
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

!!$    case('NICAM-NETCDF')
!!$
!!$       omdlid = iNICAM
!!$       if ( do_read_ocean ) call ParentOceanSetupNICAM( odims, timelen, & ! (out)
!!$                                                        basename_ocean  ) ! (in)
!!$       update_coord = .false.
!!$
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
!    case( iSCALE, iWRFARW, iNICAM )
    case( iSCALE, iWRFARW )
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
       basename_land,     &
       basename_ocean,    &
       mdlid_land, &
       mdlid_ocean,            &
       ldims,             &
       odims,             &
       use_file_landwater, &
       init_landwater_ratio, &
       init_ocean_alb_lw, &
       init_ocean_alb_sw, &
       init_ocean_z0w, &
       intrp_iter_max, &
       soilwater_ds2vc_flag, &
       elevation_collection, &
       boundary_flag,        &
       timelen,          &
       skiplen )
    use scale_comm, only: &
         COMM_bcast, &
         COMM_vars8, &
         COMM_wait
    use scale_const, only: &
         EPS => CONST_EPS, &
         UNDEF => CONST_UNDEF, &
         I_SW => CONST_I_SW, &
         I_LW => CONST_I_LW
    use scale_interp, only: &
         INTRP_factor2d, &
         INTRP_interp2d
    use scale_land_grid_cartesC, only: &
         LCZ => LAND_GRID_CARTESC_CZ
    use scale_atmos_thermodyn, only: &
         THERMODYN_specific_heat  => ATMOS_THERMODYN_specific_heat, &
         THERMODYN_rhot2temp_pres => ATMOS_THERMODYN_rhot2temp_pres
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
!!$    use mod_realinput_nicam, only: &
!!$         ParentOceanOpenNICAM, &
!!$         ParentOceanInputNICAM, &
!!$         ParentLandInputNICAM
    use mod_realinput_grads, only: &
         ParentOceanOpenGrADS, &
         ParentOceanInputGrADS, &
         ParentLandInputGrADS
    use mod_atmos_vars, only: &
         DENS, &
         MOMZ, &
         MOMX, &
         MOMY, &
         RHOT, &
         QTRC
    implicit none

    real(RP),         intent(out) :: tg  (:,:,:,:)
    real(RP),         intent(out) :: strg(:,:,:,:)
    real(RP),         intent(out) :: lst (:,:,:)
    real(RP),         intent(out) :: albg(:,:,:,:)
    real(RP),         intent(inout) :: tc_urb(IA,JA)
    real(RP),         intent(inout) :: qc_urb(IA,JA)
    real(RP),         intent(inout) :: uc_urb(IA,JA)
    real(RP),         intent(inout) :: ust   (IA,JA)
    real(RP),         intent(inout) :: albu  (IA,JA,2)
    real(RP),         intent(out) :: tw  (:,:,:)
    real(RP),         intent(out) :: sst (:,:,:)
    real(RP),         intent(out) :: albw(:,:,:,:)
    real(RP),         intent(out) :: z0w (:,:,:)
    character(len=*), intent(in)  :: basename_land
    character(len=*), intent(in)  :: basename_ocean
    integer,          intent(in)  :: mdlid_land
    integer,          intent(in)  :: mdlid_ocean
    integer,          intent(in)  :: ldims(3)
    integer,          intent(in)  :: odims(2)
    logical,          intent(in)  :: use_file_landwater   ! use land water data from files
    real(RP),         intent(in)  :: init_landwater_ratio ! Ratio of land water to storage is constant,
                                                          ! if use_file_landwater is ".false."
    real(RP),         intent(in)  :: init_ocean_alb_lw
    real(RP),         intent(in)  :: init_ocean_alb_sw
    real(RP),         intent(in)  :: init_ocean_z0w
    integer,          intent(in)  :: intrp_iter_max
    logical,          intent(in)  :: soilwater_ds2vc_flag
    logical,          intent(in)  :: elevation_collection
    logical,          intent(in)  :: boundary_flag    ! switch for making boundary file
    integer,          intent(in)  :: timelen          ! time steps in one file
    integer,          intent(in)  :: skiplen          ! skip steps

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
    real(RP) :: z0w_org  (        odims(1),odims(2))
    real(RP) :: omask    (        odims(1),odims(2))
    real(RP) :: lst_ocean(        odims(1),odims(2))

    real(RP) :: hfact_o(odims(1),odims(2),itp_nh)
    integer  :: igrd_o (odims(1),odims(2),itp_nh)
    integer  :: jgrd_o (odims(1),odims(2),itp_nh)

    real(RP) :: Qdry, Rtot, CVtot, CPtot
    real(RP) :: temp, pres

    integer :: i, j
    integer :: n, nn
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[RealinputSurface]/Categ[Input]'

    first = .true.

    if ( first ) then ! read data only once

       ! urban data

       do j = 1, JA
       do i = 1, IA
          call THERMODYN_specific_heat( QA, &
                                        qtrc(KS,i,j,:), &
                                        TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! [IN]
                                        Qdry, Rtot, CVtot, CPtot                                 ) ! [OUT]
          call THERMODYN_rhot2temp_pres( dens(KS,i,j), rhot(KS,i,j), Rtot, CVtot, CPtot, &
                                         temp, pres                                      )

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

!!$       case( iNICAM ) ! TYPE: NICAM-NETCDF
!!$
!!$          call ParentOceanOpenNICAM( olon_org, olat_org, & ! (out)
!!$                                     omask_org,          & ! (out)
!!$                                     basename_ocean,     & ! (in)
!!$                                     odims               ) ! (in)
!!$
       case( iGrADS ) ! TYPE: GrADS format

          call ParentOceanOpenGrADS

       end select

    end if


    do n = skiplen+1, timelen
       nn = n - skiplen

       if ( do_read_land ) then

          call PROF_rapstart('___SurfaceInput',3)

          select case( mdlid_land )
          case( iSCALE ) ! TYPE: SCALE-RM

             call ParentLandInputSCALE( &
                  tg_org, strg_org,           & ! (out)
                  lst_org, ust_org, albg_org, & ! (out)
                  topo_org, lmask_org,        & ! (out)
                  llon_org, llat_org, lz_org, & ! (out)
                  basename_land, ldims,       & ! (in)
                  use_file_landwater, n       ) ! (in)

          case( iWRFARW ) ! TYPE: WRF-ARW

             call ParentLandInputWRFARW( &
                  tg_org, strg_org,           & ! (out)
                  lst_org, ust_org, albg_org, & ! (out)
                  topo_org, lmask_org,        & ! (out)
                  llon_org, llat_org, lz_org, & ! (out)
                  basename_land, ldims,       & ! (in)
                  use_file_landwater, n       ) ! (in)

!!$          case( iNICAM ) ! TYPE: NICAM-NETCDF
!!$
!!$             call ParentLandInputNICAM( &
!!$                  tg_org, strg_org,           & ! (out)
!!$                  lst_org,                    & ! (out)
!!$                  llon_org, llat_org, lz_org, & ! (out)
!!$                  topo_org, lmask_org,        & ! (out)
!!$                  basename_land, ldims,       & ! (in)
!!$                  use_file_landwater, n       ) ! (in)
!!$             ust_org = UNDEF
!!$             albg_org = UNDEF
!!$
          case( iGrADS ) ! TYPE: GrADS format

             call ParentLandInputGrADS( &
                  tg_org, strg_org, smds_org, & ! (out)
                  lst_org,                    & ! (out)
                  llon_org, llat_org, lz_org, & ! (out)
                  topo_org, lmask_org,        & ! (out)
                  basename_land, ldims,       & ! (in)
                  use_file_landwater, n       ) ! (in)
             ust_org = UNDEF
             albg_org = UNDEF

          end select

          call PROF_rapend  ('___SurfaceInput',3)

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
          call COMM_bcast( lz_org, ldims(1) )
       end if

       if ( do_read_ocean ) then

          call PROF_rapstart('___SurfaceInput',3)

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

!!$          case( iNICAM ) ! TYPE: NICAM-NETCDF
!!$
!!$             call ParentOceanInputNICAM( &
!!$                  tw_org, sst_org,       & ! (out)
!!$                  basename_ocean, odims, & ! (in)
!!$                  omask_org,             & ! (in)
!!$                  n                      ) ! (in)
!!$             albw_org = UNDEF
!!$             z0w_org = UNDEF
!!$
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

          call PROF_rapend  ('___SurfaceInput',3)

       end if

       if ( serial_ocean ) then
          call PROF_rapstart('___SurfaceBcast',3)
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
          call PROF_rapend  ('___SurfaceBcast',3)
       end if

       call PROF_rapstart('___SurfaceInterp',3)

       if ( first .or. update_coord ) then
          ! interpolation factor between outer ocean grid and land grid
          call INTRP_factor2d( itp_nh,             & ! [IN]
                               ldims(2), ldims(3), & ! [IN]
                               llon_org(:,:),      & ! [IN]
                               llat_org(:,:),      & ! [IN]
                               odims(1), odims(2), & ! [IN]
                               olon_org(:,:),      & ! [IN]
                               olat_org(:,:),      & ! [IN]
                               igrd_o  (:,:,:),    & ! [OUT]
                               jgrd_o  (:,:,:),    & ! [OUT]
                               hfact_o (:,:,:)     ) ! [OUT]
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
               tg(:,:,:,nn), strg(:,:,:,nn), & ! (out)
               lst(:,:,nn), albg(:,:,:,nn),  & ! (out)
               ust, albu,                    & ! (out)
               tg_org, strg_org, smds_org,   & ! (inout)
               lst_org, albg_org,            & ! (inout)
               ust_org,                      & ! (inout)
               sst_org,                      & ! (in)
               lmask_org,                    & ! (in)
               lsmask_nest,                  & ! (in)
               topo_org,                     & ! (in)
               lz_org, llon_org, llat_org,   & ! (in)
               LCZ, LON, LAT,                & ! (in)
               ldims, odims,                 & ! (in)
               maskval_tg, maskval_strg,     & ! (in)
               init_landwater_ratio,         & ! (in)
               use_file_landwater,           & ! (in)
               use_waterratio,               & ! (in)
               soilwater_ds2vc_flag,         & ! (in)
               elevation_collection,         & ! (in)
               intrp_iter_max                ) ! (in)

       end if ! first

       if ( first .or. update_coord ) then
          ! land surface temperature at ocean grid
          call INTRP_interp2d( itp_nh,             & ! [IN]
                               ldims(2), ldims(3), & ! [IN]
                               odims(1), odims(2), & ! [IN]
                               igrd_o   (:,:,:),   & ! [IN]
                               jgrd_o   (:,:,:),   & ! [IN]
                               hfact_o  (:,:,:),   & ! [IN]
                               lst_org  (:,:),     & ! [IN]
                               lst_ocean(:,:)      ) ! [OUT]
       end if

       call replace_misval_map( sst_org, lst_ocean, odims(1), odims(2), "SST")
       call replace_misval_map( tw_org,  lst_ocean, odims(1), odims(2), "OCEAN_TEMP")

       do j = 1, odims(2)
       do i = 1, odims(1)
          if ( albw_org(i,j,I_LW) == UNDEF ) albw_org(i,j,I_LW) = init_ocean_alb_lw
          if ( albw_org(i,j,I_SW) == UNDEF ) albw_org(i,j,I_SW) = init_ocean_alb_sw
          if ( z0w_org(i,j) == UNDEF ) z0w_org(i,j) = init_ocean_z0w
       end do
       end do


       if ( first .or. update_coord ) then
          ! interporation for ocean variables
          call INTRP_factor2d( itp_nh,             & ! [IN]
                               odims(1), odims(2), & ! [IN]
                               olon_org(:,:),      & ! [IN]
                               olat_org(:,:),      & ! [IN]
                               IA, JA,             & ! [IN]
                               lon     (:,:),      & ! [IN]
                               lat     (:,:),      & ! [IN]
                               igrd    (:,:,:),    & ! [OUT]
                               jgrd    (:,:,:),    & ! [OUT]
                               hfact   (:,:,:)     ) ! [OUT]
       end if

       call INTRP_interp2d( itp_nh,               & ! [IN]
                            odims(1), odims(2),   & ! [IN]
                            IA, JA,               & ! [IN]
                            igrd    (:,:,:),      & ! [IN]
                            jgrd    (:,:,:),      & ! [IN]
                            hfact   (:,:,:),      & ! [IN]
                            tw_org  (:,:),        & ! [IN]
                            tw      (:,:,nn)      ) ! [OUT]

       call INTRP_interp2d( itp_nh,               & ! [IN]
                            odims(1), odims(2),   & ! [IN]
                            IA, JA,               & ! [IN]
                            igrd    (:,:,:),      & ! [IN]
                            jgrd    (:,:,:),      & ! [IN]
                            hfact   (:,:,:),      & ! [IN]
                            sst_org (:,:),        & ! [IN]
                            sst     (:,:,nn)      ) ! [OUT]

       call INTRP_interp2d( itp_nh,               & ! [IN]
                            odims(1), odims(2),   & ! [IN]
                            IA, JA,               & ! [IN]
                            igrd    (:,:,:),      & ! [IN]
                            jgrd    (:,:,:),      & ! [IN]
                            hfact   (:,:,:),      & ! [IN]
                            albw_org(:,:,I_LW),   & ! [IN]
                            albw    (:,:,I_LW,nn) ) ! [OUT]

       call INTRP_interp2d( itp_nh,               & ! [IN]
                            odims(1), odims(2),   & ! [IN]
                            IA, JA,               & ! [IN]
                            igrd    (:,:,:),      & ! [IN]
                            jgrd    (:,:,:),      & ! [IN]
                            hfact   (:,:,:),      & ! [IN]
                            albw_org(:,:,I_SW),   & ! [IN]
                            albw    (:,:,I_SW,nn) ) ! [OUT]

       call INTRP_interp2d( itp_nh,               & ! [IN]
                            odims(1), odims(2),   & ! [IN]
                            IA, JA,               & ! [IN]
                            igrd    (:,:,:),      & ! [IN]
                            jgrd    (:,:,:),      & ! [IN]
                            hfact   (:,:,:),      & ! [IN]
                            z0w_org (:,:),        & ! [IN]
                            z0w     (:,:,nn)      ) ! [OUT]

       if ( first ) then

          ! replace values over the ocean ####
          do j = 1, JA
          do i = 1, IA
             if( abs(lsmask_nest(i,j)-0.0_RP) < EPS ) then ! ocean grid
                lst(i,j,nn) = sst(i,j,nn)
                ust(i,j)    = sst(i,j,nn)
             endif
          enddo
          enddo

       end if

       first = .false.

       call PROF_rapend  ('___SurfaceInterp',3)

       ! required one-step data only
       if( .NOT. boundary_flag ) exit

    end do ! time loop

    return
  end subroutine ParentSurfaceInput

  !> Boundary Data Write
  subroutine ParentSurfaceBoundary( &
       tg, &
       strg, &
       lst, &
       albg, &
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
    use scale_file_cartesC, only: &
         FILE_CARTESC_create, &
         FILE_CARTESC_def_var, &
         FILE_CARTESC_enddef, &
         FILE_CARTESC_write_var
    use scale_time, only: &
         TIME_NOWDATE
    implicit none

    real(RP),         intent(in)   :: tg(:,:,:,:)
    real(RP),         intent(in)   :: strg(:,:,:,:)
    real(RP),         intent(in)   :: lst(:,:,:)
    real(RP),         intent(in)   :: albg(:,:,:,:)
    real(RP),         intent(in)   :: tw(:,:,:,:)
    real(RP),         intent(in)   :: sst(:,:,:)
    real(RP),         intent(in)   :: albw(:,:,:,:)
    real(RP),         intent(in)   :: z0(:,:,:)
    real(DP),         intent(in)   :: update_dt
    character(len=*), intent(in)   :: basename
    character(len=*), intent(in)   :: title
    integer,          intent(in)   :: numsteps ! total time steps

    character(len=H_SHORT) :: boundary_out_dtype = 'DEFAULT'  !< REAL4 or REAL8
    integer :: nowdate(6)
    integer :: fid, vid(10)
    integer :: ts, te
    !---------------------------------------------------------------------------

    call PROF_rapstart('___SurfaceOutput',3)

    ts = 1
    te = numsteps

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[RealinputSurface]/Categ[Boundary]'

    nowdate = TIME_NOWDATE
    nowdate(1) = nowdate(1)

    call FILE_CARTESC_create( basename, title, boundary_out_dtype, fid, date=nowdate )

    call FILE_CARTESC_def_var( fid,                      & ! [IN]
         'LAND_TEMP', 'Reference Land Temperature', 'K', & ! [IN]
         'LXYT', boundary_out_dtype,                     & ! [IN]
         vid(1),                                         & ! [OUT]
         timeintv=update_dt, nsteps=numsteps             ) ! [IN]
    call FILE_CARTESC_def_var( fid,                        & ! [IN]
         'LAND_WATER', 'Reference Land Moisture', 'm3/m3', & ! [IN]
         'LXYT', boundary_out_dtype,                       & ! [IN]
         vid(2),                                           & ! [OUT]
         timeintv=update_dt, nsteps=numsteps               ) ! [IN]
    call FILE_CARTESC_def_var( fid,                                  & ! [IN]
         'LAND_SFC_TEMP', 'Reference Land Surface Temperature', 'K', & ! [IN]
          'XYT', boundary_out_dtype,                                 & ! [IN]
         vid(3),                                                     & ! [OUT]
         timeintv=update_dt, nsteps=numsteps                         ) ! [IN]
    call FILE_CARTESC_def_var( fid,                                     & ! [IN]
         'LAND_ALB_LW', 'Reference Land Surface Albedo Long-wave', '1', & ! [IN]
          'XYT', boundary_out_dtype,                                    & ! [IN]
         vid(4),                                                        & ! [OUT]
         timeintv=update_dt, nsteps=numsteps                            ) ! [IN]
    call FILE_CARTESC_def_var( fid,                                      & ! [IN]
         'LAND_ALB_SW', 'Reference Land Surface Albedo Short-wave', '1', & ! [IN]
          'XYT', boundary_out_dtype,                                     & ! [IN]
         vid(5),                                                         & ! [OUT]
         timeintv=update_dt, nsteps=numsteps                             ) ! [IN]
    call FILE_CARTESC_def_var( fid,                        & ! [IN]
         'OCEAN_TEMP', 'Reference Ocean Temperature', 'K', & ! [IN]
          'OXYT', boundary_out_dtype,                      & ! [IN]
         vid(6),                                           & ! [OUT]
         timeintv=update_dt, nsteps=numsteps               ) ! [IN]
    call FILE_CARTESC_def_var( fid,                                    & ! [IN]
         'OCEAN_SFC_TEMP', 'Reference Ocean Surface Temperature', 'K', & ! [IN]
          'XYT', boundary_out_dtype,                                   & ! [IN]
         vid(7),                                                       & ! [OUT]
         timeintv=update_dt, nsteps=numsteps                           ) ! [IN]
    call FILE_CARTESC_def_var( fid,                                       & ! [IN]
         'OCEAN_ALB_LW', 'Reference Ocean Surface Albedo Long-wave', '1', & ! [IN]
          'XYT', boundary_out_dtype,                                      & ! [IN]
         vid(8),                                                          & ! [OUT]
         timeintv=update_dt, nsteps=numsteps                              ) ! [IN]
    call FILE_CARTESC_def_var( fid,                                        & ! [IN]
         'OCEAN_ALB_SW', 'Reference Ocean Surface Albedo Short-wave', '1', & ! [IN]
          'XYT', boundary_out_dtype,                                       & ! [IN]
         vid(9),                                                           & ! [OUT]
         timeintv=update_dt, nsteps=numsteps                               ) ! [IN]
    call FILE_CARTESC_def_var( fid,                         & ! [IN]
         'OCEAN_SFC_Z0', 'Reference Ocean Surface Z0', 'm', & ! [IN]
          'XYT', boundary_out_dtype,                        & ! [IN]
         vid(10),                                           & ! [OUT]
         timeintv=update_dt, nsteps=numsteps                ) ! [IN]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var( fid, vid(1),  tg  (:,:,:,     ts:te), 'LAND_TEMP',      'LXYT', update_dt )
    call FILE_CARTESC_write_var( fid, vid(2),  strg(:,:,:,     ts:te), 'LAND_WATER',     'LXYT', update_dt )
    call FILE_CARTESC_write_var( fid, vid(3),  lst (  :,:,     ts:te), 'LAND_SFC_TEMP',  'XYT',  update_dt )
    call FILE_CARTESC_write_var( fid, vid(4),  albg(  :,:,I_LW,ts:te), 'LAND_ALB_LW',    'XYT',  update_dt )
    call FILE_CARTESC_write_var( fid, vid(5),  albg(  :,:,I_SW,ts:te), 'LAND_ALB_SW',    'XYT',  update_dt )
    call FILE_CARTESC_write_var( fid, vid(6),  tw  (:,:,:,     ts:te), 'OCEAN_TEMP',     'OXYT', update_dt )
    call FILE_CARTESC_write_var( fid, vid(7),  sst (  :,:,     ts:te), 'OCEAN_SFC_TEMP', 'XYT',  update_dt )
    call FILE_CARTESC_write_var( fid, vid(8),  albw(  :,:,I_LW,ts:te), 'OCEAN_ALB_LW',   'XYT',  update_dt )
    call FILE_CARTESC_write_var( fid, vid(9),  albw(  :,:,I_SW,ts:te), 'OCEAN_ALB_SW',   'XYT',  update_dt )
    call FILE_CARTESC_write_var( fid, vid(10), z0  (  :,:,     ts:te), 'OCEAN_SFC_Z0',   'XYT',  update_dt )

    call PROF_rapend  ('___SurfaceOutput',3)

    return
  end subroutine ParentSurfaceBoundary


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
    use scale_interp, only: &
         INTRP_factor2d, &
         INTRP_factor3d, &
         INTRP_interp2d, &
         INTRP_interp3d
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
    call INTRP_factor2d( itp_nh,             & ! [IN]
                         odims(1), odims(2), & ! [IN]
                         olon_org(:,:),      & ! [IN]
                         olat_org(:,:),      & ! [IN]
                         ldims(2), ldims(3), & ! [IN]
                         llon_org(:,:),      & ! [IN]
                         llat_org(:,:),      & ! [IN]
                         igrd_l  (:,:,:),    & ! [OUT]
                         jgrd_l  (:,:,:),    & ! [OUT]
                         hfact_l (:,:,:)     ) ! [OUT]

    ! sst on land grid
    call INTRP_interp2d( itp_nh,             & ! [IN]
                         odims(1), odims(2), & ! [IN]
                         ldims(2), ldims(3), & ! [IN]
                         igrd_l   (:,:,:),   & ! [IN]
                         jgrd_l   (:,:,:),   & ! [IN]
                         hfact_l  (:,:,:),   & ! [IN]
                         sst_org  (:,:),     & ! [IN]
                         sst_land (:,:)      ) ! [OUT]

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

    call INTRP_factor3d( itp_nh,                & ! [IN]
                         ldims(1), 1, ldims(1), & ! [IN]
                         ldims(2), ldims(3),    & ! [IN]
                         llon_org(:,:),         & ! [IN]
                         llat_org(:,:),         & ! [IN]
                         lz3d_org(:,:,:),       & ! [IN]
                         LKMAX, LKS, LKE,       & ! [IN]
                         IA, JA,                & ! [IN]
                         lon     (:,:),         & ! [IN]
                         lat     (:,:),         & ! [IN]
                         lcz_3D  (:,:,:),       & ! [IN]
                         igrd    (    :,:,:),   & ! [OUT]
                         jgrd    (    :,:,:),   & ! [OUT]
                         hfact   (    :,:,:),   & ! [OUT]
                         kgrdl   (:,:,:,:,:),   & ! [OUT]
                         vfactl  (:,:,:,:,:)    ) ! [OUT]

    call INTRP_interp2d( itp_nh,             & ! [IN]
                         ldims(2), ldims(3), & ! [IN]
                         IA, JA,             & ! [IN]
                         igrd    (:,:,:),    & ! [IN]
                         jgrd    (:,:,:),    & ! [IN]
                         hfact   (:,:,:),    & ! [IN]
                         lst_org (:,:),      & ! [IN]
                         lst     (:,:)       ) ! [OUT]

    call INTRP_interp2d( itp_nh,             & ! [IN]
                         ldims(2), ldims(3), & ! [IN]
                         IA, JA,             & ! [IN]
                         igrd    (:,:,:),    & ! [IN]
                         jgrd    (:,:,:),    & ! [IN]
                         hfact   (:,:,:),    & ! [IN]
                         ust_org (:,:),      & ! [IN]
                         ust     (:,:)       ) ! [OUT]

    call INTRP_interp2d( itp_nh,             & ! [IN]
                         ldims(2), ldims(3), & ! [IN]
                         IA, JA,             & ! [IN]
                         igrd    (:,:,:),    & ! [IN]
                         jgrd    (:,:,:),    & ! [IN]
                         hfact   (:,:,:),    & ! [IN]
                         albg_org(:,:,I_LW), & ! [IN]
                         albg    (:,:,I_LW)  ) ! [OUT]

    call INTRP_interp2d( itp_nh,             & ! [IN]
                         ldims(2), ldims(3), & ! [IN]
                         IA, JA,             & ! [IN]
                         igrd    (:,:,:),    & ! [IN]
                         jgrd    (:,:,:),    & ! [IN]
                         hfact   (:,:,:),    & ! [IN]
                         albg_org(:,:,I_SW), & ! [IN]
                         albg    (:,:,I_SW)  ) ! [OUT]

    call INTRP_interp3d( itp_nh,                       & ! [IN]
                         ldims(1), ldims(2), ldims(3), & ! [IN]
                         LKMAX, LKS, LKE,              & ! [IN]
                         IA, JA,                       & ! [IN]
                         igrd  (    :,:,:),            & ! [IN]
                         jgrd  (    :,:,:),            & ! [IN]
                         hfact (    :,:,:),            & ! [IN]
                         kgrdl (:,:,:,:,:),            & ! [IN]
                         vfactl(:,:,:,:,:),            & ! [IN]
                         tg_org(:,:,:),                & ! [IN]
                         tg    (:,:,:)                 ) ! [OUT]

    do j = 1, JA
    do i = 1, IA
       tg(LKMAX,i,j) = tg(LKMAX-1,i,j)
    enddo
    enddo

    ! replace values over the ocean
    do k = 1, LKMAX
       call replace_misval_const( tg(k,:,:), maskval_tg, lsmask_nest )
    enddo


    ! elevation collection
    if ( elevation_collection ) then
       call INTRP_interp2d( itp_nh,             & ! [IN]
                            ldims(2), ldims(3), & ! [IN]
                            IA, JA,             & ! [IN]
                            igrd    (:,:,:),    & ! [IN]
                            jgrd    (:,:,:),    & ! [IN]
                            hfact   (:,:,:),    & ! [IN]
                            topo_org(:,:),      & ! [IN]
                            topo    (:,:)       ) ! [OUT]

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

          call INTRP_interp3d( itp_nh,                       & ! [IN]
                               ldims(1), ldims(2), ldims(3), & ! [IN]
                               LKMAX, LKS, LKE,              & ! [IN]
                               IA, JA,                       & ! [IN]
                               igrd    (    :,:,:),          & ! [IN]
                               jgrd    (    :,:,:),          & ! [IN]
                               hfact   (    :,:,:),          & ! [IN]
                               kgrdl   (:,:,:,:,:),          & ! [IN]
                               vfactl  (:,:,:,:,:),          & ! [IN]
                               smds_org(:,:,:),              & ! [IN]
                               smds    (:,:,:)               ) ! [OUT]

          do k = 1, LKMAX-1
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

          call INTRP_interp3d( itp_nh,                       & ! [IN]
                               ldims(1), ldims(2), ldims(3), & ! [IN]
                               LKMAX, LKS, LKE,              & ! [IN]
                               IA, JA,                       & ! [IN]
                               igrd    (    :,:,:),          & ! [IN]
                               jgrd    (    :,:,:),          & ! [IN]
                               hfact   (    :,:,:),          & ! [IN]
                               kgrdl   (:,:,:,:,:),          & ! [IN]
                               vfactl  (:,:,:,:,:),          & ! [IN]
                               strg_org(:,:,:),              & ! [IN]
                               strg    (:,:,:)               ) ! [OUT]
       end if

       ! replace values over the ocean
       do k = 1, LKMAX-1
          call replace_misval_const( strg(k,:,:), maskval_strg, lsmask_nest )
       enddo

       do j = 1, JA
       do i = 1, IA
          strg(LKMAX,i,j) = strg(LKMAX-1,i,j)
       enddo
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
      data,     &
      lsmask,   &
      nx,       &
      ny,       &
      landdata, &
      iter_max  )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    integer,  intent(in)    :: nx
    integer,  intent(in)    :: ny
    real(RP), intent(inout) :: data  (nx,ny)
    real(RP), intent(in)    :: lsmask(nx,ny)
    logical,  intent(in)    :: landdata   ! .true. => land data , .false. => ocean data
    integer,  intent(in)    :: iter_max

    integer  :: mask     (nx,ny)
    integer  :: mask_prev(nx,ny)
    real(RP) :: data_prev(nx,ny)
    real(RP) :: tmp, cnt, sw
    integer  :: mask_target

    integer  :: num_land, num_ocean, num_replaced
    integer  :: istr, iend, jstr, jend
    integer  :: i, j, ii, jj, ite
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [interp_OceanLand_data]/Categ[realinit]'

    if ( landdata ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** target mask : LAND'
       mask_target = 1 ! interpolation for land data
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** target mask : OCEAN'
       mask_target = 0 ! interpolation for ocean data
    endif

    ! search target cell for interpolation
    num_land  = 0
    num_ocean = 0
    do j = 1, ny
    do i = 1, nx
       mask(i,j) = int( 0.5_RP - sign(0.5_RP,abs(lsmask(i,j)-1.0_RP)-EPS) ) ! 1 for land, 0 for ocean
       num_land  = num_land  + (   mask(i,j) )
       num_ocean = num_ocean + ( 1-mask(i,j) )
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,3I8,A,I8)') '*** ite = ', 0, &
               ', (land,ocean,replaced) = ', num_land, num_ocean, 0, ' / ', nx*ny

    ! start interpolation
    do ite = 1, iter_max
       ! save previous state
       mask_prev(:,:) = mask(:,:)
       data_prev(:,:) = data(:,:)
       num_replaced   = 0

       do j  = 1, ny
       do i  = 1, nx

          if( mask(i,j) == mask_target ) cycle ! already filled

          ! collect neighbor grid
          istr = max(i-1,1 )
          iend = min(i+1,nx)
          jstr = max(j-1,1 )
          jend = min(j+1,ny)

          tmp = 0.0_RP
          cnt = 0.0_RP
          do jj = jstr, jend
          do ii = istr, iend
             sw = 0.5_RP - sign(0.5_RP,real(abs(mask_prev(ii,jj)-mask_target),kind=RP)-EPS)

             tmp = tmp + sw * data_prev(ii,jj)
             cnt = cnt + sw
          enddo
          enddo

          if ( cnt >= 3.0_RP ) then ! replace by average of neighbor grid value
             data(i,j) = tmp / cnt
             mask(i,j) = mask_target

             num_replaced = num_replaced + 1
          endif

       enddo
       enddo

       if ( landdata ) then
          num_land  = num_land  + num_replaced
          num_ocean = num_ocean - num_replaced
       else
          num_land  = num_land  - num_replaced
          num_ocean = num_ocean + num_replaced
       endif
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,3I8,A,I8)') '*** ite = ', ite, &
                  ', (land,ocean,replaced) = ', num_land, num_ocean, num_replaced, ' / ', nx*ny

       if( num_replaced == 0 ) exit

    enddo ! itelation


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
             write(*,*) "xxx data for mask of "//trim(elem)//"(",i,",",j,") includes missing value."
             write(*,*) "xxx Please check input data of SKINTEMP or SST. "
             call PRC_MPIstop
          else
             data(i,j) = maskval(i,j)
          endif
       endif
    enddo
    enddo

  end subroutine replace_misval_map

end module mod_realinput
