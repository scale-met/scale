!-------------------------------------------------------------------------------
!> module REAL input
!!
!! @par Description
!!          read data from file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module mod_realinput
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_ocean_grid_cartesC_index
  use scale_land_grid_cartesC_index
  use scale_urban_grid_cartesC_index
  use scale_index
  use scale_tracer
  use scale_cpl_sfc_index

  use scale_prc, only: &
     PRC_IsMaster, &
     PRC_abort
  use scale_comm_cartesC, only: &
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
  private :: ParentAtmosFinalize
  private :: ParentAtmosOpen
  private :: ParentAtmosInput
  private :: BoundaryAtmosSetup
  private :: BoundaryAtmosOutput

  private :: ParentSurfaceSetup
  private :: ParentSurfaceFinalize
  private :: ParentSurfaceOpen
  private :: ParentSurfaceInput
  private :: BoundarySurfaceSetup
  private :: BoundarySurfaceOutput
  private :: interp_OceanLand_data

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, public, parameter :: iNetCDF  = 1
  integer, public, parameter :: iGrADS  = 2

  integer, private :: KA_org, KS_org, KE_org
  integer, private :: IA_org, IS_org, IE_org
  integer, private :: JA_org, JS_org, JE_org
  integer, private :: lKA_org, lKS_org, lKE_org
  integer, private :: lIA_org, lIS_org, lIE_org
  integer, private :: lJA_org, lJS_org, lJE_org
  integer, private :: oIA_org, oIS_org, oIE_org
  integer, private :: oJA_org, oJS_org, oJE_org

  real(RP), private, allocatable :: LON_org (:,:)
  real(RP), private, allocatable :: LAT_org (:,:)
  real(RP), private, allocatable :: CZ_org  (:,:,:)

  real(RP), private, allocatable :: W_org   (:,:,:) ! scalar point
  real(RP), private, allocatable :: U_org   (:,:,:) ! scalar point
  real(RP), private, allocatable :: V_org   (:,:,:) ! scalar point
  real(RP), private, allocatable :: DENS_org(:,:,:)
  real(RP), private, allocatable :: PT_org  (:,:,:)
  real(RP), private, allocatable :: TEMP_org(:,:,:)
  real(RP), private, allocatable :: PRES_org(:,:,:)
  real(RP), private, allocatable :: QTRC_org (:,:,:,:)
  real(RP), private, allocatable :: QV_org   (:,:,:)
  real(RP), private, allocatable :: RH_org   (:,:,:)
  real(RP), private, allocatable :: QHYD_org (:,:,:,:)
  real(RP), private, allocatable :: QNUM_org (:,:,:,:)

  real(RP), private, allocatable :: llon_org (:,:)
  real(RP), private, allocatable :: llat_org (:,:)

  real(RP), private, allocatable :: olon_org (:,:)
  real(RP), private, allocatable :: olat_org (:,:)

  real(RP), private, allocatable :: rotc_org(:,:,:)

  integer,  private              :: itp_nh_a  = 4 ! for atmos
  integer,  private              :: itp_nh_l  = 4 ! for land
  integer,  private              :: itp_nh_o  = 4 ! for ocean
  integer,  private              :: itp_nh_ol = 5 ! for ocean-land

  integer,  private, parameter   :: I_intrp_linear = 0
  integer,  private, parameter   :: I_intrp_dstwgt = 1
  integer,  private              :: itp_type_a
  integer,  private              :: itp_type_l
  integer,  private              :: itp_type_o

  integer,  private, allocatable :: igrd (    :,:,:)
  integer,  private, allocatable :: jgrd (    :,:,:)
  real(RP), private, allocatable :: hfact(    :,:,:)
  integer,  private, allocatable :: kgrd (:,:,:,:,:)
  real(RP), private, allocatable :: vfact(:,  :,:,:)

  integer,  private, allocatable :: oigrd (:,:,:)
  integer,  private, allocatable :: ojgrd (:,:,:)
  real(RP), private, allocatable :: ohfact(:,:,:)

  logical,  private              :: ol_interp
  real(RP), private, allocatable :: hfact_ol(:,:,:)
  integer,  private, allocatable :: igrd_ol (:,:,:)
  integer,  private, allocatable :: jgrd_ol (:,:,:)

  logical,  private              :: serial_atmos
  logical,  private              :: serial_land
  logical,  private              :: serial_ocean
  logical,  private              :: do_read_atmos
  logical,  private              :: do_read_land
  logical,  private              :: do_read_ocean

  logical,  private              :: mixing_ratio
  logical,  private              :: temp2pt
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

  integer,                private :: FILTER_ORDER               = 8       ! order of the hyper-diffusion (must be even)
  integer,                private :: FILTER_NITER               = 0       ! times for hyper-diffusion iteration

  logical,                private :: USE_FILE_DENSITY           = .false. ! use density data from files
  logical,                private :: pt_dry                     = .true.  ! potential temperature is defined associated with dry air
  logical,                private :: SAME_MP_TYPE               = .false. ! microphysics type of the parent model is same as it in this model
  character(len=4),       private :: upper_qv_type = "ZERO" !< how qv is given at higher level than outer model
                                                            !< "ZERO": 0
                                                            !< "COPY": copy values from the highest level of outer model

  character(len=H_SHORT), private :: INTRP_TYPE                 = "LINEAR" ! "LINEAR" or "DIST-WEIGHT"
                                                                           !   LINEAR     : bi-linear interpolation
                                                                           !   DIST-WEIGHT: distance-weighted mean of the nearest N-neighbors

  logical, private :: first_atmos   = .true.
  logical, private :: first_surface = .true.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine REALINPUT_atmos
    use scale_const, only: &
       P00 => CONST_PRE00
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
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_specific_heat
    implicit none

    logical :: USE_SFC_DIAGNOSES          = .false. !> use surface diagnoses
    logical :: USE_DATA_UNDER_SFC         = .true.  !> use data under the surface
    logical :: USE_NONHYDRO_DENS_BOUNDARY = .false. !> use non-hydrostatic density for boundary data
    logical :: SKIP_VERTICAL_RANGE_CHECK  = .false. !> skip chkecking if the domain top does not exceed that of the parent in the vertical direction


    namelist / PARAM_MKINIT_REAL_ATMOS / &
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
       FILTER_ORDER,               &
       FILTER_NITER,               &
       USE_FILE_DENSITY,           &
       PT_dry,                     &
       USE_NONHYDRO_DENS_BOUNDARY, &
       USE_SFC_DIAGNOSES,          &
       USE_DATA_UNDER_SFC,         &
       SAME_MP_TYPE,               &
       UPPER_QV_TYPE,              &
       INTRP_TYPE,                 &
       SKIP_VERTICAL_RANGE_CHECK

    character(len=6) :: basename_num
    character(len=H_LONG) :: basename_out_mod
    character(len=19)     :: timelabel

    integer  :: dims(6) ! dims 1-3: normal, 4-6: staggerd
    integer  :: timelen

    integer  :: fid_atmos
    integer  :: vid_atmos(5+QA)
    logical  :: qtrc_flag(QA)

    real(RP) :: DENS_in(KA,IA,JA)
    real(RP) :: MOMZ_in(KA,IA,JA) ! staggered point
    real(RP) :: MOMX_in(KA,IA,JA) ! staggered point
    real(RP) :: MOMY_in(KA,IA,JA) ! staggered point
    real(RP) :: RHOT_in(KA,IA,JA)
    real(RP) :: QTRC_in(KA,IA,JA,QA)

    real(RP) :: VELZ_in(KA,IA,JA) ! staggered point
    real(RP) :: VELX_in(KA,IA,JA) ! staggered point
    real(RP) :: VELY_in(KA,IA,JA) ! staggered point
    real(RP) :: PT_in  (KA,IA,JA)
    real(RP) :: PRES_in(KA,IA,JA)

    real(RP) :: Qdry (KA,IA,JA)
    real(RP) :: Rtot (KA,IA,JA)
    real(RP) :: CPtot(KA,IA,JA)
    real(RP) :: CVtot(KA,IA,JA)

    integer  :: ifile, istep, t, tall
    integer  :: k, i, j, iq
    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO('REALINPUT_atmos',*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_ATMOS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("REALINPUT_atmos",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("REALINPUT_atmos",*) 'Not appropriate names in namelist PARAM_MKINIT_REAL_ATMOS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_REAL_ATMOS)

    if ( BOUNDARY_UPDATE_DT <= 0.0_DP ) then
       LOG_ERROR("REALINPUT_atmos",*) 'BOUNDARY_UPDATE_DT is necessary in real case preprocess'
       call PRC_abort
    endif

    if ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM ) then
       basename_num = '_00000'
    else
       basename_num = ''
    endif

    select case( INTRP_TYPE )
    case ( "LINEAR" )
       itp_nh_a = 4
       itp_type_a = I_intrp_linear
    case ( "DIST-WEIGHT" )
       itp_nh_a = COMM_CARTESC_NEST_INTERP_LEVEL
       itp_type_a = I_intrp_dstwgt
    case default
       LOG_ERROR("REALINPUT_atmos",*) 'Unsupported type of INTRP_TYPE : ', trim(INTRP_TYPE)
       LOG_ERROR_CONT(*) '       It must be "LINEAR" or "DIST-WEIGHT"'
       call PRC_abort
    end select

    call ParentAtmosSetup( FILETYPE_ORG,     & ! [IN]
                           basename_org,     & ! [IN]
                           basename_num,     & ! [IN]
                           SERIAL_PROC_READ, & ! [IN]
                           SAME_MP_TYPE,     & ! [IN]
                           USE_FILE_DENSITY, & ! [IN]
                           dims(:),          & ! [OUT]
                           timelen,          & ! [OUT]
                           qtrc_flag(:)      ) ! [OUT]

    if ( timelen > 0 ) then
       NUMBER_OF_TSTEPS = timelen ! read from file
    endif

    LOG_NEWLINE
    LOG_INFO("REALINPUT_atmos",*) 'Number of temporal data in each file : ', NUMBER_OF_TSTEPS

    do ifile = 1, NUMBER_OF_FILES

       if ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM ) then
          write(basename_num,'(A,I5.5)') '_', ifile-1
       else
          basename_num = ''
       endif

       LOG_NEWLINE
       LOG_INFO("REALINPUT_atmos",*) 'read external data from : ', trim(basename_org)//trim(basename_num)

       call ParentAtmosOpen( FILETYPE_ORG, basename_org, basename_num ) ![IN]

       do istep = 1, NUMBER_OF_TSTEPS

          tall = NUMBER_OF_TSTEPS * (ifile-1) + istep ! consecutive time step (input)
          t    = tall - NUMBER_OF_SKIP_TSTEPS         ! time step (output)

          if ( t <= 0 ) then
             LOG_PROGRESS('(1x,A,I4,A,I5,A,I6,A)') &
                          '[file,step,cons.] = [', ifile, ',', istep, ',', tall, '] ...skip.'
             cycle
          endif

          if ( t == 1 .OR. BASENAME_BOUNDARY /= '' ) then

             LOG_PROGRESS('(1x,A,I4,A,I5,A,I6,A)') &
                          '[file,step,cons.] = [', ifile, ',', istep, ',', tall, ']'

             ! read prepared data
             call ParentAtmosInput( FILETYPE_ORG,              & ! [IN]
                                    basename_num,              & ! [IN]
                                    dims(:),                   & ! [IN]
                                    istep,                     & ! [IN]
                                    USE_SFC_DIAGNOSES,         & ! [IN]
                                    USE_DATA_UNDER_SFC,        & ! [IN]
                                    SAME_MP_TYPE,              & ! [IN]
                                    SKIP_VERTICAL_RANGE_CHECK, & ! [IN]
                                    qtrc_flag(:),              & ! [IN]
                                    DENS_in(:,:,:),            & ! [OUT]
                                    MOMZ_in(:,:,:),            & ! [OUT]
                                    MOMX_in(:,:,:),            & ! [OUT]
                                    MOMY_in(:,:,:),            & ! [OUT]
                                    RHOT_in(:,:,:),            & ! [OUT]
                                    QTRC_in(:,:,:,:),          & ! [OUT]
                                    VELZ_in(:,:,:),            & ! [OUT]
                                    VELX_in(:,:,:),            & ! [OUT]
                                    VELY_in(:,:,:),            & ! [OUT]
                                    PT_in  (:,:,:),            & ! [OUT]
                                    PRES_in(:,:,:)             ) ! [OUT]
          else
             LOG_PROGRESS('(1x,A,I4,A,I5,A,I6,A)') &
                          '[file,step,cons.] = [', ifile, ',', istep, ',', tall, '] ...skip.'
          endif

          !--- store prognostic variables as initial
          if ( t == 1 ) then
             LOG_NEWLINE
             LOG_INFO("REALINPUT_atmos",*) 'store initial state.'

             !$omp parallel do collapse(2)
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

             !$omp parallel do collapse(3)
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
                                         qtrc_flag(:),       & ! [IN]
                                         fid_atmos,          & ! [OUT]
                                         vid_atmos(:)        ) ! [OUT]
             endif

             if ( use_nonhydro_dens_boundary ) then
                call ATMOS_THERMODYN_specific_heat( KA, KS, KE, IA, 1, IA, JA, 1, JA, QA, &
                                                    QTRC_in(:,:,:,:),                                        & ! [IN]
                                                    TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! [IN]
                                                    Qdry(:,:,:), Rtot(:,:,:), CVtot(:,:,:), CPtot(:,:,:)     ) ! [OUT]
                !$omp parallel do collapse(2)
                do j = 1, JA
                do i = 1, IA
                do k = KS, KE
                   DENS_in(k,i,j) = ( PRES_in(k,i,j) / P00 )**( CVtot(k,i,j) / CPtot(k,i,j) ) * P00 / ( Rtot(k,i,j) * PT_in(k,i,j) )
                end do
                end do
                end do
             end if

             call BoundaryAtmosOutput( DENS_in(:,:,:),     & ! [IN]
                                       VELZ_in(:,:,:),     & ! [IN]
                                       VELX_in(:,:,:),     & ! [IN]
                                       VELY_in(:,:,:),     & ! [IN]
                                       PT_in  (:,:,:),     & ! [IN]
                                       QTRC_in(:,:,:,:),   & ! [IN]
                                       qtrc_flag(:),       & ! [IN]
                                       fid_atmos,          & ! [IN]
                                       vid_atmos(:),       & ! [IN]
                                       BOUNDARY_UPDATE_DT, & ! [IN]
                                       t                   ) ! [IN]
          endif

       enddo ! istep loop
    enddo ! ifile loop

    call ParentAtmosFinalize( FILETYPE_ORG ) ![IN]

    return
  end subroutine REALINPUT_atmos

  !-----------------------------------------------------------------------------
  subroutine REALINPUT_surface
    use scale_const, only: &
       TEM00 => CONST_TEM00
    use scale_time, only: &
       TIME_gettimelabel
    use scale_landuse, only: &
       fact_ocean => LANDUSE_fact_ocean, &
       fact_land  => LANDUSE_fact_land,  &
       fact_urban => LANDUSE_fact_urban, &
       LANDUSE_PFT_nmax
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_SFC_TEMP,   &
       ATMOS_PHY_SF_SFC_albedo, &
       ATMOS_PHY_SF_SFC_Z0M,    &
       ATMOS_PHY_SF_SFC_Z0H,    &
       ATMOS_PHY_SF_SFC_Z0E
    use mod_ocean_admin, only: &
       OCEAN_do
    use scale_ocean_phy_ice_simple, only: &
       OCEAN_PHY_ICE_freezetemp
    use mod_ocean_vars, only: &
       ICE_FLAG,         &
       OCEAN_TEMP,       &
       OCEAN_SALT,       &
       OCEAN_UVEL,       &
       OCEAN_VVEL,       &
       OCEAN_OCN_Z0M,    &
       OCEAN_ICE_TEMP,   &
       OCEAN_ICE_MASS,   &
       OCEAN_SFC_TEMP,   &
       OCEAN_SFC_albedo, &
       OCEAN_SFC_Z0M,    &
       OCEAN_SFC_Z0H,    &
       OCEAN_SFC_Z0E
    use mod_land_admin, only: &
       LAND_do
    use mod_land_vars, only: &
       LAND_TEMP,       &
       LAND_WATER,      &
       LAND_ICE,        &
       LAND_SFC_TEMP,   &
       LAND_SFC_albedo
    use mod_urban_admin, only: &
       URBAN_do
    use mod_urban_vars, only: &
       URBAN_TC,        &
       URBAN_QC,        &
       URBAN_UC,        &
       URBAN_TR,        &
       URBAN_TB,        &
       URBAN_TG,        &
       URBAN_TRL,       &
       URBAN_TBL,       &
       URBAN_TGL,       &
       URBAN_RAINR,     &
       URBAN_RAINB,     &
       URBAN_RAING,     &
       URBAN_ROFF,      &
       URBAN_SFC_TEMP,  &
       URBAN_SFC_albedo
    implicit none

    logical                  :: USE_FILE_LANDWATER                          = .true.    ! use land water data from files
    real(RP)                 :: INIT_LANDWATER_RATIO                        = 0.5_RP    ! Ratio of land water to storage is constant, if USE_FILE_LANDWATER is ".false." (all PFT)
!    real(RP)                 :: INIT_LANDWATER_RATIO_EACH(LANDUSE_PFT_nmax)             ! Ratio of land water to storage is constant, if USE_FILE_LANDWATER is ".false." (each PFT)
    real(RP)                 :: INIT_OCEAN_ALB_LW                           = 0.04_RP   ! initial LW albedo on the ocean
    real(RP)                 :: INIT_OCEAN_ALB_SW                           = 0.10_RP   ! initial SW albedo on the ocean
    real(RP)                 :: INIT_OCEAN_Z0W                              = 1.0E-3_RP ! initial surface roughness on the ocean
    character(len=H_SHORT)   :: INTRP_LAND_TEMP                             = 'off'
    character(len=H_SHORT)   :: INTRP_LAND_WATER                            = 'off'
    character(len=H_SHORT)   :: INTRP_LAND_SFC_TEMP                         = 'off'
    character(len=H_SHORT)   :: INTRP_OCEAN_TEMP                            = 'off'
    character(len=H_SHORT)   :: INTRP_OCEAN_SFC_TEMP                        = 'off'
    integer                  :: INTRP_ITER_MAX                              = 5
    character(len=H_SHORT)   :: SOILWATER_DS2VC                             = 'limit'
    logical                  :: soilwater_DS2VC_flag                                  ! true: 'critical', false: 'limit'
    logical                  :: elevation_correction                        = .true.
    logical                  :: elevation_correction_land
    logical                  :: elevation_correction_ocean

    namelist / PARAM_MKINIT_REAL_LAND / &
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
!       INIT_LANDWATER_RATIO_EACH,  &
       INTRP_TYPE,                 &
       INTRP_LAND_TEMP,            &
       INTRP_LAND_WATER,           &
       INTRP_LAND_SFC_TEMP,        &
       INTRP_ITER_MAX,             &
       FILTER_ORDER,               &
       FILTER_NITER,               &
       SOILWATER_DS2VC,            &
       ELEVATION_CORRECTION,       &
       SERIAL_PROC_READ

    namelist / PARAM_MKINIT_REAL_OCEAN / &
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
       INTRP_TYPE,                 &
       INTRP_OCEAN_TEMP,           &
       INTRP_OCEAN_SFC_TEMP,       &
       INTRP_ITER_MAX,             &
       FILTER_ORDER,               &
       FILTER_NITER,               &
       SERIAL_PROC_READ

    character(len=6) :: FILETYPE_LAND
    character(len=6) :: FILETYPE_OCEAN
    character(len=6) :: BASENAME_LAND
    character(len=6) :: BASENAME_OCEAN

    ! land
    real(RP) :: LAND_TEMP_in      (LKMAX,IA,JA)
    real(RP) :: LAND_WATER_in     (LKMAX,IA,JA)
    real(RP) :: LAND_SFC_TEMP_in  (      IA,JA)
    real(RP) :: LAND_SFC_albedo_in(      IA,JA,N_RAD_DIR,N_RAD_RGN)

    ! urban
    real(RP) :: URBAN_TC_in(IA,JA)
    real(RP) :: URBAN_QC_in(IA,JA)
    real(RP) :: URBAN_UC_in(IA,JA)
    real(RP) :: URBAN_SFC_TEMP_in(IA,JA)
    real(RP) :: URBAN_SFC_albedo_in(IA,JA,N_RAD_DIR,N_RAD_RGN)

    ! ocean
    real(RP) :: OCEAN_TEMP_in      (OKMAX,IA,JA)
    real(RP) :: OCEAN_SFC_TEMP_in  (      IA,JA)
    real(RP) :: OCEAN_SFC_albedo_in(      IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP) :: OCEAN_SFC_Z0_in    (      IA,JA)

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
    logical               :: BASENAME_ADD_NUM_LAND
    logical               :: BASENAME_ADD_NUM_OCEAN

    real(RP), allocatable :: lst_ocean(:,:)
    real(RP), allocatable :: lz_org(:)
    real(RP), allocatable :: topo_org(:,:)
    real(RP), allocatable :: lmask_org(:,:)
    real(RP), allocatable :: omask_org(:,:)

    integer :: mdlid_land, mdlid_ocean
    integer :: ldims(3), odims(2)

    integer :: totaltimesteps = 1
    integer :: timelen
    integer :: ierr

    character(len=H_LONG) :: basename_out_mod
    character(len=19)     :: timelabel

    integer  :: fid_sfc
    integer  :: vid_sfc(10)

    logical :: land_flag
    logical :: multi_land
    logical :: multi_ocean

    integer :: oistep, listep
    integer :: tall, t
    integer :: idir, irgn

    integer :: k, i, j, n
    !---------------------------------------------------------------------------

    if ( LAND_do .or. URBAN_do ) then
       land_flag = .true.
    else
       land_flag = .false.
    end if

    if ( .not. land_flag .or. .not. OCEAN_do ) then
       LOG_ERROR("REALINPUT_surface",*) 'OCEAN_ and LAND_DYN_TYPE must be set'
    end if


    LOG_NEWLINE
    LOG_INFO('REALINPUT_surface',*) 'Setup LAND'

    ! LAND/URBAN
!    INIT_LANDWATER_RATIO_EACH(:) = -1.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_LAND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("REALINPUT_surface",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("REALINPUT_surface",*) 'Not appropriate names in namelist PARAM_MKINIT_REAL_LAND. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_REAL_LAND)

    NUMBER_OF_FILES_LAND            = NUMBER_OF_FILES
    NUMBER_OF_TSTEPS_LAND           = NUMBER_OF_TSTEPS
    NUMBER_OF_SKIP_TSTEPS_LAND      = NUMBER_OF_SKIP_TSTEPS
    FILETYPE_LAND                   = FILETYPE_ORG
    BASENAME_ADD_NUM_LAND           = BASENAME_ADD_NUM
    BASENAME_BOUNDARY_LAND          = BASENAME_BOUNDARY
    BOUNDARY_POSTFIX_TIMELABEL_LAND = BOUNDARY_POSTFIX_TIMELABEL
    BOUNDARY_TITLE_LAND             = BOUNDARY_TITLE
    BOUNDARY_UPDATE_DT_LAND         = BOUNDARY_UPDATE_DT
    elevation_correction_land       = elevation_correction

    if ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM_LAND ) then
       BASENAME_LAND = "_00000"
    else
       BASENAME_LAND = ""
    endif

!!$    if( .NOT. USE_FILE_LANDWATER ) then
!!$       if( all( INIT_LANDWATER_RATIO_EACH(:) < 0.0_RP ) ) then
!!$          LOG_INFO("REALINPUT_surface",*) 'Applied INIT_LANDWATER_RATIO, instead of INIT_LANDWATER_RATIO_EACH.'
!!$          INIT_LANDWATER_RATIO_EACH(:) = INIT_LANDWATER_RATIO
!!$       else
!!$          if( any( INIT_LANDWATER_RATIO_EACH(:) < 0.0_RP ) ) then
!!$             LOG_ERROR("REALINPUT_surface",*) 'Insufficient elemtents of array (INIT_LANDWATER_RATIO_EACH):', INIT_LANDWATER_RATIO_EACH(:)
!!$             call PRC_abort
!!$          endif
!!$       endif
!!$    endif

    select case( SOILWATER_DS2VC )
    case( 'critical' )
       SOILWATER_DS2VC_flag = .true.
    case('limit' )
       SOILWATER_DS2VC_flag = .false.
    case default
       LOG_ERROR("REALINPUT_surface",*) 'Unsupported SOILWATER_DS2CV TYPE:', trim(SOILWATER_DS2VC)
       call PRC_abort
    end select

    !serial_land = SERIAL_PROC_READ
    ! if serial_land = .false., interp_OceanLand_data may result in different value
    serial_land = .true.

    select case( INTRP_TYPE )
    case ( "LINEAR" )
       itp_nh_l  = 4
       itp_type_l = I_intrp_linear
    case ( "DIST-WEIGHT" )
       itp_nh_l  = COMM_CARTESC_NEST_INTERP_LEVEL
       itp_type_l = I_intrp_dstwgt
    case default
       LOG_ERROR("REALINPUT_surface",*) 'Unsupported type of INTRP_TYPE : ', trim(INTRP_TYPE)
       LOG_ERROR_CONT(*) '       It must be "LINEAR" or "DIST-WEIGHT"'
       call PRC_abort
    end select



    LOG_NEWLINE
    LOG_INFO('REALINPUT_surface',*) 'Setup OCEAN'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_OCEAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("REALINPUT_surface",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("REALINPUT_surface",*) 'Not appropriate names in namelist PARAM_MKINIT_REAL_OCEAN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_REAL_OCEAN)

    NUMBER_OF_FILES_OCEAN            = NUMBER_OF_FILES
    NUMBER_OF_TSTEPS_OCEAN           = NUMBER_OF_TSTEPS
    NUMBER_OF_SKIP_TSTEPS_OCEAN      = NUMBER_OF_SKIP_TSTEPS
    FILETYPE_OCEAN                   = FILETYPE_ORG
    BASENAME_ADD_NUM_OCEAN           = BASENAME_ADD_NUM
    BASENAME_BOUNDARY_OCEAN          = BASENAME_BOUNDARY
    BOUNDARY_POSTFIX_TIMELABEL_OCEAN = BOUNDARY_POSTFIX_TIMELABEL
    BOUNDARY_TITLE_OCEAN             = BOUNDARY_TITLE
    BOUNDARY_UPDATE_DT_OCEAN         = BOUNDARY_UPDATE_DT
    elevation_correction_ocean       = elevation_correction

    if ( NUMBER_OF_FILES > 1 .OR. BASENAME_ADD_NUM_OCEAN ) then
       BASENAME_OCEAN = "_00000"
    else
       BASENAME_OCEAN = ""
    endif

    ! if serial_ocean = .false., interp_OceanLand_data may result in different value
    !serial_ocean = SERIAL_PROC_READ
    serial_ocean = .true.

    select case( INTRP_TYPE )
    case ( "LINEAR" )
       itp_nh_o  = 4
       itp_type_o = I_intrp_linear
    case ( "DIST-WEIGHT" )
       itp_nh_o  = COMM_CARTESC_NEST_INTERP_LEVEL
       itp_type_o = I_intrp_dstwgt
    case default
       LOG_ERROR("REALINPUT_surface",*) 'Unsupported type of INTRP_TYPE : ', trim(INTRP_TYPE)
       LOG_ERROR_CONT(*) '       It must be "LINEAR" or "DIST-WEIGHT"'
       call PRC_abort
    end select

    itp_nh_ol = COMM_CARTESC_NEST_INTERP_LEVEL

    multi_land  = ( NUMBER_OF_FILES_LAND * NUMBER_OF_TSTEPS_LAND - NUMBER_OF_SKIP_TSTEPS_LAND ) > 1
    multi_ocean = BASENAME_BOUNDARY_OCEAN .ne. ''

    if ( ( multi_land .and. multi_ocean ) .AND. &
       ( ( NUMBER_OF_FILES_LAND            .NE.   NUMBER_OF_FILES_OCEAN            ) .OR. &
         ( NUMBER_OF_TSTEPS_LAND           .NE.   NUMBER_OF_TSTEPS_OCEAN           ) .OR. &
         ( NUMBER_OF_SKIP_TSTEPS_LAND      .NE.   NUMBER_OF_SKIP_TSTEPS_OCEAN      ) .OR. &
         ( BASENAME_BOUNDARY_LAND          .NE.   BASENAME_BOUNDARY_OCEAN          ) .OR. &
         ( BOUNDARY_POSTFIX_TIMELABEL_LAND .NEQV. BOUNDARY_POSTFIX_TIMELABEL_OCEAN ) .OR. &
         ( BOUNDARY_TITLE_LAND             .NE.   BOUNDARY_TITLE_OCEAN             ) .OR. &
         ( BOUNDARY_UPDATE_DT_LAND         .NE.   BOUNDARY_UPDATE_DT_OCEAN         ) ) ) then
       LOG_ERROR("REALINPUT_surface",*) 'The following LAND/OCEAN parameters must be consistent due to technical problem:'
       LOG_ERROR_CONT(*) '           NUMBER_OF_FILES, NUMBER_OF_TSTEPS, NUMBER_OF_SKIP_TSTEPS,'
       LOG_ERROR_CONT(*) '           BASENAME_BOUNDARY, BOUNDARY_POSTFIX_TIMELABEL, BOUNDARY_TITLE, BOUNDARY_UPDATE_DT.'
       call PRC_abort
    end if

    if ( multi_land .and. .not. multi_ocean ) then
       LOG_ERROR("REALINPUT_surface",*) 'To output boundary data for land data, it is necessary to output ocean data'
       call PRC_abort
    end if

    call ParentSurfaceSetup( ldims, odims,           & ![OUT]
                             mdlid_land,             & ![OUT]
                             mdlid_ocean,            & ![OUT]
                             timelen,                & ![OUT]
                             BASENAME_ORG,           & ![IN]
                             BASENAME_LAND,          & ![IN]
                             BASENAME_OCEAN,         & ![IN]
                             FILETYPE_LAND,          & ![IN]
                             FILETYPE_OCEAN,         & ![IN]
                             USE_FILE_LANDWATER,     & ![IN]
                             INTRP_ITER_MAX,         & ![IN]
                             intrp_land_TEMP,        & ![IN]
                             intrp_land_water,       & ![IN]
                             intrp_land_sfc_TEMP,    & ![IN]
                             intrp_ocean_TEMP,       & ![IN]
                             intrp_ocean_sfc_TEMP    ) ![IN]

    if ( timelen > 0 ) then
       NUMBER_OF_TSTEPS_OCEAN = timelen ! read from file
    endif

    totaltimesteps = NUMBER_OF_FILES_OCEAN * NUMBER_OF_TSTEPS_OCEAN

    if ( totaltimesteps <= NUMBER_OF_SKIP_TSTEPS_OCEAN ) then
       LOG_ERROR("REALSINPUT_surface",*) "NUMBER_OF_FILES_OCEAN * NUMBER_OF_TSTEPS_OCEAN <= NUMBER_OF_SKIP_TSTEPS_OCEAN"
       call PRC_abort
    end if
    if ( multi_ocean ) then
       if ( BOUNDARY_UPDATE_DT_OCEAN <= 0.0_DP ) then
          LOG_ERROR("REALINPUT_surface",*) 'BOUNDARY_UPDATE_DT is necessary in real case preprocess'
          call PRC_abort
       endif
    end if

    allocate( lz_org(lKA_org) )
    allocate( topo_org(lIA_org, lJA_org) )
    allocate( lmask_org(lIA_org, lJA_org) )
    allocate( omask_org(oIA_org, oJA_org) )
    allocate( lst_ocean(oIA_org, oJA_org) )

    !--- read external file
    do n = 1, NUMBER_OF_FILES_OCEAN

       if ( NUMBER_OF_FILES_LAND > 1 .OR. BASENAME_ADD_NUM_LAND ) then
          write(BASENAME_LAND,'(A,I5.5)') '_', n-1
       else
          BASENAME_LAND = ''
       endif
       if ( NUMBER_OF_FILES_OCEAN > 1 .OR. BASENAME_ADD_NUM_OCEAN ) then
          write(BASENAME_OCEAN,'(A,I5.5)') '_', n-1
       else
          BASENAME_OCEAN = ''
       endif

       LOG_NEWLINE
       LOG_INFO("REALINPUT_surface",*) 'Target Postfix Name (Land) : ', trim(BASENAME_LAND)
       LOG_INFO("REALINPUT_surface",*) 'Target Postfix Name (Ocean): ', trim(BASENAME_OCEAN)
       LOG_INFO("REALINPUT_surface",*) 'Time Steps in One File     : ', NUMBER_OF_TSTEPS

       call ParentSurfaceOpen( FILETYPE_LAND, FILETYPE_OCEAN, & ! [IN]
                               basename_org, basename_land, basename_ocean ) ! [IN]

       do oistep = 1, NUMBER_OF_TSTEPS_OCEAN

          tall = NUMBER_OF_TSTEPS_OCEAN * (n-1) + oistep
          t    = tall - NUMBER_OF_SKIP_TSTEPS_OCEAN
          if ( t <= 0 ) then
             LOG_PROGRESS('(1x,A,I4,A,I5,A,I6,A)') &
                  '[file,step,cons.] = [', n, ',', oistep, ',', tall, '] ...skip.'
             cycle
          end if
          if ( t > 1 .and. .not. multi_ocean .and. .not. multi_land ) exit

          if ( multi_land ) then
             listep = oistep
          else
             listep = NUMBER_OF_SKIP_TSTEPS_LAND + 1
          end if

          ! read all prepared data
          call ParentSurfaceInput( LAND_TEMP_in       (:,:,:),           &
                                   LAND_WATER_in      (:,:,:),           &
                                   LAND_SFC_TEMP_in   (:,:),             &
                                   LAND_SFC_albedo_in (:,:,:,:),         &
                                   URBAN_TC_in(:,:),                     &
                                   URBAN_QC_in(:,:),                     &
                                   URBAN_UC_in(:,:),                     &
                                   URBAN_SFC_TEMP_in(:,:),               &
                                   URBAN_SFC_albedo_in(:,:,:,:),         &
                                   lst_ocean(:,:),                       &
                                   lz_org(:), topo_org(:,:),             &
                                   lmask_org(:,:), omask_org(:,:),       &
                                   OCEAN_TEMP_in      (OKS,:,:),         &
                                   OCEAN_SFC_TEMP_in  (    :,:),         &
                                   OCEAN_SFC_albedo_in(    :,:,:,:),     &
                                   OCEAN_SFC_Z0_in    (    :,:),         &
                                   BASENAME_ORG,                         &
                                   BASENAME_LAND, BASENAME_OCEAN,        &
                                   mdlid_land, mdlid_ocean,              &
                                   ldims, odims,                         &
                                   USE_FILE_LANDWATER,                   &
                                   INIT_LANDWATER_RATIO,                 &
!                                   INIT_LANDWATER_RATIO_EACH(:),         &
                                   INIT_OCEAN_ALB_LW, INIT_OCEAN_ALB_SW, &
                                   INIT_OCEAN_Z0W,                       &
                                   INTRP_ITER_MAX,                       &
                                   SOILWATER_DS2VC_flag,                 &
                                   elevation_correction_land,            &
                                   elevation_correction_ocean,           &
                                   oistep, listep,                       &
                                   multi_land,                           &
                                   URBAN_do                              )

          !--- input initial data
          if ( t == 1 ) then

             !$omp parallel do
             do j = 1, JA
             do i = 1, IA
                OCEAN_SFC_TEMP(i,j) = OCEAN_SFC_TEMP_in(i,j)
                OCEAN_SFC_Z0M (i,j) = OCEAN_SFC_Z0_in  (i,j)
                OCEAN_SFC_Z0H (i,j) = OCEAN_SFC_Z0_in  (i,j)
                OCEAN_SFC_Z0E (i,j) = OCEAN_SFC_Z0_in  (i,j)
                do irgn = I_R_IR, I_R_VIS
                   do idir = I_R_direct, I_R_diffuse
                      OCEAN_SFC_albedo(i,j,idir,irgn) = OCEAN_SFC_albedo_in(i,j,idir,irgn)
                   enddo
                enddo
                do k = 1, OKMAX
                   OCEAN_TEMP(k,i,j) = OCEAN_TEMP_in(OKS,i,j)
                   OCEAN_SALT(k,i,j) = 0.0_RP
                   OCEAN_UVEL(k,i,j) = 0.0_RP
                   OCEAN_VVEL(k,i,j) = 0.0_RP
                enddo
                OCEAN_OCN_Z0M (i,j) = OCEAN_SFC_Z0_in(i,j)
                if ( ICE_FLAG ) then
                   OCEAN_ICE_TEMP(i,j) = min( OCEAN_SFC_TEMP_in(i,j), OCEAN_PHY_ICE_freezetemp )
                   OCEAN_ICE_MASS(i,j) = 0.0_RP
                end if

                LAND_SFC_TEMP(i,j) = LAND_SFC_TEMP_in(i,j)
                do irgn = I_R_IR, I_R_VIS
                do idir = I_R_direct, I_R_diffuse
                   LAND_SFC_albedo(i,j,idir,irgn) = LAND_SFC_albedo_in(i,j,idir,irgn)
                enddo
                enddo
                do k = 1, LKMAX
                   LAND_TEMP (k,i,j) = LAND_TEMP_in(k,i,j)
                   if ( LAND_TEMP(k,i,j) >= TEM00 ) then
                      LAND_WATER(k,i,j) = LAND_WATER_in(k,i,j)
                      LAND_ICE  (k,i,j) = 0.0_RP
                   else
                      LAND_WATER(k,i,j) = 0.0_RP
                      LAND_ICE(k,i,j)   = LAND_WATER_in(k,i,j)
                   end if
                enddo

                if ( URBAN_do ) then
                   URBAN_SFC_TEMP(i,j) = URBAN_SFC_TEMP_in(i,j)
                   do irgn = I_R_IR, I_R_VIS
                   do idir = I_R_direct, I_R_diffuse
                      URBAN_SFC_albedo(i,j,idir,irgn) = URBAN_SFC_albedo_in(i,j,idir,irgn)
                   enddo
                   enddo
                   do k = 1, UKMAX
                      URBAN_TRL(k,i,j) = URBAN_SFC_TEMP_in(i,j)
                      URBAN_TBL(k,i,j) = URBAN_SFC_TEMP_in(i,j)
                      URBAN_TGL(k,i,j) = URBAN_SFC_TEMP_in(i,j)
                   enddo
                   URBAN_TC   (i,j) = URBAN_TC_in      (i,j)
                   URBAN_QC   (i,j) = URBAN_QC_in      (i,j)
                   URBAN_UC   (i,j) = URBAN_UC_in      (i,j)
                   URBAN_TR   (i,j) = URBAN_SFC_TEMP_in(i,j)
                   URBAN_TB   (i,j) = URBAN_SFC_TEMP_in(i,j)
                   URBAN_TG   (i,j) = URBAN_SFC_TEMP_in(i,j)
                   URBAN_RAINR(i,j) = 0.0_RP
                   URBAN_RAINB(i,j) = 0.0_RP
                   URBAN_RAING(i,j) = 0.0_RP
                   URBAN_ROFF (i,j) = 0.0_RP
                end if

                ATMOS_PHY_SF_SFC_Z0M (i,j) = OCEAN_SFC_Z0M(i,j)
                ATMOS_PHY_SF_SFC_Z0H (i,j) = OCEAN_SFC_Z0H(i,j)
                ATMOS_PHY_SF_SFC_Z0E (i,j) = OCEAN_SFC_Z0E(i,j)

                if ( URBAN_do ) then
                   ATMOS_PHY_SF_SFC_TEMP(i,j) = fact_ocean(i,j) * OCEAN_SFC_TEMP(i,j) &
                                              + fact_land (i,j) * LAND_SFC_TEMP (i,j) &
                                     + fact_urban(i,j) * URBAN_SFC_TEMP(i,j)
                   do irgn = I_R_IR, I_R_VIS
                   do idir = I_R_direct, I_R_diffuse
                      ATMOS_PHY_SF_SFC_albedo(i,j,idir,irgn) = fact_ocean(i,j) * OCEAN_SFC_albedo(i,j,idir,irgn) &
                                                             + fact_land (i,j) * LAND_SFC_albedo (i,j,idir,irgn) &
                                                             + fact_urban(i,j) * URBAN_SFC_albedo(i,j,idir,irgn)
                   enddo
                   enddo
                else
                   ATMOS_PHY_SF_SFC_TEMP(i,j) = fact_ocean(i,j) * OCEAN_SFC_TEMP(i,j) &
                                              + fact_land (i,j) * LAND_SFC_TEMP (i,j)
                   do irgn = I_R_IR, I_R_VIS
                   do idir = I_R_direct, I_R_diffuse
                      ATMOS_PHY_SF_SFC_albedo(i,j,idir,irgn) = fact_ocean(i,j) * OCEAN_SFC_albedo(i,j,idir,irgn) &
                                                             + fact_land (i,j) * LAND_SFC_albedo (i,j,idir,irgn)
                   enddo
                   enddo
                endif
             enddo
             enddo

          end if ! t==1


          !--- output boundary data
          if ( multi_ocean ) then

             if ( t==1 ) then
                if ( BOUNDARY_POSTFIX_TIMELABEL_OCEAN ) then
                   call TIME_gettimelabel( timelabel )
                   basename_out_mod = trim(BASENAME_BOUNDARY_OCEAN)//'_'//trim(timelabel)
                else
                   basename_out_mod = trim(BASENAME_BOUNDARY_OCEAN)
                endif

                call BoundarySurfaceSetup( basename_out_mod,         & ! [IN]
                                           BOUNDARY_TITLE_OCEAN,     & ! [IN]
                                           BOUNDARY_UPDATE_DT_OCEAN, & ! [IN]
                                           multi_ocean, multi_land,  & ! [IN]
                                           fid_sfc,                  & ! [OUT]
                                           vid_sfc(:)                ) ! [OUT]
             end if ! t==1

             call BoundarySurfaceOutput( LAND_TEMP_in(:,:,:),      & ! [IN]
                                         LAND_WATER_in(:,:,:),     & ! [IN]
                                         LAND_SFC_TEMP_in(:,:),    & ! [IN]
                                         OCEAN_TEMP_in(:,:,:),     & ! [IN]
                                         OCEAN_SFC_TEMP_in(:,:),   & ! [IN]
                                         OCEAN_SFC_Z0_in(:,:),     & ! [IN]
                                         multi_ocean, multi_land,  & ! [IN]
                                         fid_sfc,                  & ! [IN]
                                         vid_sfc(:),               & ! [IN]
                                         BOUNDARY_UPDATE_DT_OCEAN, & ! [IN]
                                         t                         ) ! [IN]
          end if

       end do ! oistep

       ! required one-step data only
       if( .not. ( multi_land .or. multi_ocean ) ) exit

    enddo

    deallocate( lst_ocean )
    deallocate( lz_org )
    deallocate( topo_org )
    deallocate( lmask_org )
    deallocate( omask_org )

    call ParentSurfaceFinalize( FILETYPE_LAND, & ![IN]
                                FILETYPE_OCEAN ) ![IN]

    return
  end subroutine REALINPUT_surface


  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtmosSetup( &
       inputtype,           &
       basename,            &
       basename_num,        &
       serial_in,           &
       same_mptype,         &
       use_file_density_in, &
       dims,                &
       timelen,             &
       qtrc_flag            )
    use scale_comm_cartesC, only: &
       COMM_world, &
       COMM_datatype
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use scale_mapprojection, only: &
       mappingparam, &
       mappinginfo, &
       MAPPROJECTION_rotcoef, &
       MAPPROJECTION_get_param
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LON, &
       ATMOS_GRID_CARTESC_REAL_LAT
    use mod_realinput_netcdf, only: &
       ParentAtmosSetupNetCDF
    use mod_realinput_grads, only: &
       ParentAtmosSetupGrADS
    implicit none

    character(len=*), intent(in)  :: inputtype
    character(len=*), intent(in)  :: basename
    character(len=*), intent(in)  :: basename_num
    logical,          intent(in)  :: serial_in           ! read by a serial process
    logical,          intent(in)  :: same_mptype
    logical,          intent(in)  :: use_file_density_in ! use density data from files

    integer, intent(out) :: dims(6)
    integer, intent(out) :: timelen
    logical, intent(out) :: qtrc_flag(QA)

    real(RP), pointer :: LON_all(:,:)
    real(RP), pointer :: LAT_all(:,:)

    real(RP) :: LON_min, LON_max
    real(RP) :: LAT_min, LAT_max

    type(mappingparam) :: mapping_param
    type(mappinginfo)  :: mapping_info

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    serial_atmos = serial_in
    if ( serial_atmos .and. (.not. PRC_IsMaster) ) then
       do_read_atmos = .false.
    else
       do_read_atmos = .true.
    endif

    select case(inputtype)
    case('NetCDF')

       call ParentAtmosSetupNetCDF( dims(:),      & ! [OUT]
                                    timelen,      & ! [OUT]
                                    mixing_ratio, & ! [OUT]
                                    update_coord, & ! [OUT]
                                    mapping_info, & ! [OUT]
                                    qtrc_flag(:), & ! [OUT]
                                    LON_all,      & ! [OUT]
                                    LAT_all,      & ! [OUT]
                                    basename,     & ! [IN]
                                    basename_num, & ! [IN]
                                    same_mptype,  & ! [IN]
                                    pt_dry,       & ! [INOUT]
                                    serial_atmos, & ! [INOUT]
                                    do_read_atmos ) ! [INOUT]

    case('GrADS')

       if ( do_read_atmos ) then
          call ParentAtmosSetupGrADS( dims(:),      & ! [OUT]
                                      timelen,      & ! [OUT]
                                      qtrc_flag(:), & ! [OUT]
                                      LON_all,      & ! [OUT]
                                      LAT_all,      & ! [OUT]
                                      basename,     & ! [IN]
                                      basename_num  ) ! [IN]
       endif
       mapping_info%mapping_name = ""
       mixing_ratio = .false.
       update_coord = .true.

    case default

       LOG_ERROR("ParentAtmosSetup",*) 'Unsupported type of input data : ', trim(inputtype)
       call PRC_abort

    end select

    LON_min = minval( ATMOS_GRID_CARTESC_REAL_LON(:,:) )
    LON_max = maxval( ATMOS_GRID_CARTESC_REAL_LON(:,:) )
    LAT_min = minval( ATMOS_GRID_CARTESC_REAL_LAT(:,:) )
    LAT_max = maxval( ATMOS_GRID_CARTESC_REAL_LAT(:,:) )

    if ( serial_atmos ) then
       call COMM_bcast( 6, dims(:) )
       call COMM_bcast( timelen )

       call COMM_bcast( QA, qtrc_flag(:) )

       call COMM_bcast( mixing_ratio )
       call COMM_bcast( update_coord )

       call COMM_bcast( mapping_info%mapping_name )
       call COMM_bcast( mapping_info%longitude_of_central_meridian )
       call COMM_bcast( mapping_info%longitude_of_projection_origin )
       call COMM_bcast( mapping_info%latitude_of_projection_origin )
       call COMM_bcast( 2, mapping_info%standard_parallel )
       call COMM_bcast( mapping_info%rotation )

       if ( .not. do_read_atmos ) then
          allocate( LON_all(dims(2), dims(3)) )
          allocate( LAT_all(dims(2), dims(3)) )
       end if
       call COMM_bcast( dims(2), dims(3), LON_all )
       call COMM_bcast( dims(2), dims(3), LAT_all )

       call MPI_AllReduce( MPI_IN_PLACE, LON_min, 1, COMM_datatype, MPI_MIN, COMM_WORLD, ierr )
       call MPI_AllReduce( MPI_IN_PLACE, LON_max, 1, COMM_datatype, MPI_MAX, COMM_WORLD, ierr )
       call MPI_ALLReduce( MPI_IN_PLACE, LAT_min, 1, COMM_datatype, MPI_MIN, COMM_WORLD, ierr )
       call MPI_ALLReduce( MPI_IN_PLACE, LAT_max, 1, COMM_datatype, MPI_MAX, COMM_WORLD, ierr )
    endif

    call get_IJrange( IS_org, IE_org, JS_org, JE_org,     & ! [OUT]
                      dims(2), dims(3),                   & ! [IN]
                      LON_min, LON_max, LAT_min, LAT_max, & ! [IN]
                      LON_all, LAT_all                    ) ! [IN]
    KA_org = dims(1) + 2
    KS_org = 1
    KE_org = KA_org
    IA_org = IE_org - IS_org + 1
    JA_org = JE_org - JS_org + 1

    allocate( LON_org (        IA_org, JA_org) )
    allocate( LAT_org (        IA_org, JA_org) )
    allocate( CZ_org  (KA_org, IA_org, JA_org) )

    allocate( W_org   (KA_org, IA_org, JA_org    ) )
    allocate( U_org   (KA_org, IA_org, JA_org    ) )
    allocate( V_org   (KA_org, IA_org, JA_org    ) )
    allocate( PT_org  (KA_org, IA_org, JA_org    ) )
    allocate( TEMP_org(KA_org, IA_org, JA_org    ) )
    allocate( PRES_org(KA_org, IA_org, JA_org    ) )
    allocate( DENS_org(KA_org, IA_org, JA_org    ) )
    allocate( QTRC_org(KA_org, IA_org, JA_org, QA) )

    allocate( QV_org  (KA_org, IA_org, JA_org       ) )
    allocate( RH_org  (KA_org, IA_org, JA_org       ) )
    allocate( QHYD_org(KA_org, IA_org, JA_org, N_HYD) )
    allocate( QNUM_org(KA_org, IA_org, JA_org, N_HYD) )

    allocate( rotc_org(IA_org, JA_org, 2) )

    allocate( igrd (       IA, JA, itp_nh_a) )
    allocate( jgrd (       IA, JA, itp_nh_a) )
    allocate( hfact(       IA, JA, itp_nh_a) )
    allocate( kgrd (KA, 2, IA, JA, itp_nh_a) )
    allocate( vfact(KA,    IA, JA, itp_nh_a) )

    !$omp parallel do
    do j = 1, JA_org
    do i = 1, IA_org
       LON_org(i,j) = LON_all(i-1+IS_org,j-1+JS_org)
       LAT_org(i,j) = LAT_all(i-1+IS_org,j-1+JS_org)
    end do
    end do

    deallocate( LON_all )
    deallocate( LAT_all )

    if ( mapping_info%mapping_name /= "" ) then
       call MAPPROJECTION_get_param( mapping_info, mapping_param )
       call MAPPROJECTION_rotcoef( &
            IA_org, 1, IA_org, JA_org, 1, JA_org, &
            lon_org(:,:), lat_org(:,:), &
            mapping_info%mapping_name, &
            mapping_param, &
            rotc_org(:,:,1), rotc_org(:,:,2) )
    else
       rotc_org(:,:,1) = 1.0_RP
       rotc_org(:,:,2) = 0.0_RP
    end if

    first_atmos = .true.

    return
  end subroutine ParentAtmosSetup

  !-----------------------------------------------------------------------------
  !> Atmos Finalize
  subroutine ParentAtmosFinalize( &
       inputtype           )
    use mod_realinput_netcdf, only: &
       ParentAtmosFinalizeNetCDF
    implicit none

    character(len=*), intent(in)  :: inputtype
    !---------------------------------------------------------------------------

    select case(inputtype)
    case('NetCDF')

       if ( do_read_atmos ) then
          call ParentAtmosFinalizeNetCDF
       end if

    case('GrADS')
       ! do nothing
    end select

    deallocate( igrd )
    deallocate( jgrd )
    deallocate( hfact )
    deallocate( kgrd  )
    deallocate( vfact )

    deallocate( LON_org )
    deallocate( LAT_org )
    deallocate( CZ_org  )

    deallocate( W_org    )
    deallocate( U_org    )
    deallocate( V_org    )
    deallocate( PT_org   )
    deallocate( TEMP_org )
    deallocate( PRES_org )
    deallocate( DENS_org )
    deallocate( QTRC_org )

    deallocate( QV_org    )
    deallocate( RH_org    )
    deallocate( QHYD_org  )
    deallocate( QNUM_org  )

    return
  end subroutine ParentAtmosFinalize


 !-----------------------------------------------------------------------------
  !> Atmosphere Data Open
  subroutine ParentAtmosOpen( &
       inputtype, &
       basename_org, &
       basename_num )
    use mod_realinput_netcdf, only: &
       ParentAtmosOpenNetCDF
    implicit none

    character(len=*), intent(in)  :: inputtype
    character(len=*), intent(in)  :: basename_org
    character(len=*), intent(in)  :: basename_num

    !---------------------------------------------------------------------------

    select case(inputtype)
    case('NetCDF')

       if ( do_read_atmos ) then
          call ParentAtmosOpenNetCDF( basename_org, basename_num )
       endif

    case('GrADS')

       ! do nothing
    end select

    return
  end subroutine ParentAtmosOpen

  !-----------------------------------------------------------------------------
  !> Atmosphere Data Read
  subroutine ParentAtmosInput( &
       inputtype,     &
       basename_num,  &
       dims,          &
       istep,         &
       sfc_diagnoses, &
       under_sfc,     &
       same_mptype,   &
       skip_vcheck,   &
       qtrc_flag,     &
       DENS,          &
       MOMZ,          &
       MOMX,          &
       MOMY,          &
       RHOT,          &
       QTRC,          &
       VELZ,          &
       VELX,          &
       VELY,          &
       POTT,          &
       PRES           )
    use scale_const, only: &
       UNDEF   => CONST_UNDEF, &
       PI      => CONST_PI, &
       GRAV    => CONST_GRAV, &
       EPSvap  => CONST_EPSvap, &
       EPSTvap => CONST_EPSTvap, &
       PRE00   => CONST_PRE00, &
       Rdry    => CONST_Rdry, &
       CPdry   => CONST_CPdry, &
       CPvap   => CONST_CPvap, &
       LAPS    => CONST_LAPS
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_grid_cartesC_metric, only: &
       ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       N_HYD, &
       I_QV, &
       QLS, &
       QLE
    use scale_atmos_thermodyn, only: &
       THERMODYN_qdry           => ATMOS_THERMODYN_qdry, &
       THERMODYN_r              => ATMOS_THERMODYN_r,    &
       THERMODYN_cp             => ATMOS_THERMODYN_cp,   &
       THERMODYN_rhot2pres      => ATMOS_THERMODYN_rhot2pres, &
       THERMODYN_temp_pres2pott => ATMOS_THERMODYN_temp_pres2pott, &
       THERMODYN_specific_heat  => ATMOS_THERMODYN_specific_heat
    use scale_atmos_saturation, only: &
       psat => ATMOS_SATURATION_psat_liq
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_interp, only: &
       INTERP_domain_compatibility, &
       INTERP_factor3d,             &
       INTERP_interp3d
    use scale_filter, only: &
       FILTER_hyperdiff
    use mod_realinput_netcdf, only: &
       ParentAtmosInputNetCDF
    use mod_realinput_grads, only: &
       ParentAtmosInputGrADS
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    use mod_atmos_phy_sf_vars, only: &
       Z0M => ATMOS_PHY_SF_SFC_Z0M
    use scale_atmos_grid_cartesC, only: &
       CX => ATMOS_GRID_CARTESC_CX, &
       CY => ATMOS_GRID_CARTESC_CY
    use scale_atmos_grid_cartesC_real, only: &
       CZ => ATMOS_GRID_CARTESC_REAL_CZ
    use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy
    use scale_topography, only: &
       topo => TOPOGRAPHY_Zsfc
    implicit none

    character(len=*), intent(in) :: inputtype
    character(len=*), intent(in) :: basename_num
    integer,          intent(in) :: dims(6)
    integer,          intent(in) :: istep
    logical,          intent(in) :: sfc_diagnoses
    logical,          intent(in) :: under_sfc
    logical,          intent(in) :: same_mptype   ! Is microphysics type same between outer and inner model
    logical,          intent(in) :: skip_vcheck
    logical,          intent(in) :: qtrc_flag(QA)

    real(RP), intent(out) :: DENS(KA,IA,JA)
    real(RP), intent(out) :: MOMZ(KA,IA,JA)
    real(RP), intent(out) :: MOMX(KA,IA,JA)
    real(RP), intent(out) :: MOMY(KA,IA,JA)
    real(RP), intent(out) :: RHOT(KA,IA,JA)
    real(RP), intent(out) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(out) :: VELZ(KA,IA,JA)
    real(RP), intent(out) :: VELX(KA,IA,JA)
    real(RP), intent(out) :: VELY(KA,IA,JA)
    real(RP), intent(out) :: POTT(KA,IA,JA)
    real(RP), intent(out) :: PRES(KA,IA,JA)

    real(RP) :: PRES2(KA,IA,JA)
    real(RP) :: TEMP (KA,IA,JA)
    real(RP) :: W    (KA,IA,JA)
    real(RP) :: U    (KA,IA,JA)
    real(RP) :: V    (KA,IA,JA)
    real(RP) :: QV   (KA,IA,JA)
    real(RP) :: QC   (KA,IA,JA)
    real(RP) :: DENS2(KA,IA,JA)
    real(RP) :: u_on_map, v_on_map

    real(RP) :: qdry(KA_org,IA_org,JA_org)
    real(RP) :: Rtot(KA_org,IA_org,JA_org)
    real(RP) :: CPtot(KA_org,IA_org,JA_org)
    real(RP) :: CVtot(KA_org,IA_org,JA_org)
    real(RP) :: qtot, p_sat, qm
    real(RP) :: rhot_tmp

    integer :: lm_layer(IA_org,JA_org)

    real(RP) :: X_org(IA_org,JA_org)
    real(RP) :: Y_org(IA_org,JA_org)
    logical  :: zonal, pole

    real(RP) :: wsum(KA,IA,JA)
    real(RP) :: work(KA,IA,JA)

    logical :: same_mptype_
    logical :: uvmet
    logical :: nopres
    logical :: nodens
    logical :: rh2qv
    logical :: qnum_flag

    real(RP) :: one(KA,IA,JA)

    real(RP) :: dz

    integer :: kref, k0
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('___AtmosInput',3)

    if ( do_read_atmos ) then
       select case(inputtype)
       case('NetCDF')
          call ParentAtmosInputNetCDF( KA_org, KS_org, KE_org, &
                                       IA_org, IS_org, IE_org, &
                                       JA_org, JS_org, JE_org, &
                                       QA, &
                                       CZ_org  (:,:,:),   & ! [INOUT]
                                       W_org   (:,:,:),   & ! [OUT]
                                       U_org   (:,:,:),   & ! [OUT]
                                       V_org   (:,:,:),   & ! [OUT]
                                       PRES_org(:,:,:),   & ! [OUT]
                                       DENS_org(:,:,:),   & ! [OUT]
                                       TEMP_org(:,:,:),   & ! [OUT]
                                       PT_org  (:,:,:),   & ! [OUT]
                                       QTRC_org(:,:,:,:), & ! [OUT]
                                       QV_org  (:,:,:),   & ! [OUT]
                                       RH_org  (:,:,:),   & ! [OUT]
                                       QHYD_org(:,:,:,:), & ! [OUT]
                                       QNUM_org(:,:,:,:), & ! [OUT]
                                       nopres, nodens,    & ! [OUT]
                                       uvmet,             & ! [OUT]
                                       temp2pt, rh2qv,    & ! [OUT]
                                       qnum_flag,         & ! [OUT]
                                       same_mptype,       & ! [IN]
                                       sfc_diagnoses,     & ! [IN]
                                       update_coord,      & ! [IN]
                                       dims(:),           & ! [IN]
                                       istep              ) ! [IN]
          same_mptype_ = same_mptype
       case('GrADS')
          call ParentAtmosInputGrADS ( KA_org, KS_org, KE_org, &
                                       IA_org, IS_org, IE_org, &
                                       JA_org, JS_org, JE_org, &
                                       QA, &
                                       W_org   (:,:,:),        & ! [OUT]
                                       U_org   (:,:,:),        & ! [OUT]
                                       V_org   (:,:,:),        & ! [OUT]
                                       PRES_org(:,:,:),        & ! [OUT]
                                       DENS_org(:,:,:),        & ! [OUT]
                                       PT_org  (:,:,:),        & ! [OUT]
                                       TEMP_org(:,:,:),        & ! [OUT]
                                       QV_org  (:,:,:),        & ! [OUT]
                                       RH_org  (:,:,:),        & ! [OUT]
                                       QHYD_org(:,:,:,:),      & ! [OUT]
                                       QTRC_org(:,:,:,:),      & ! [OUT]
                                       CZ_org  (:,:,:),        & ! [OUT]
                                       nopres, nodens,         & ! [OUT]
                                       temp2pt, rh2qv,         & ! [OUT]
                                       basename_num,           & ! [IN]
                                       sfc_diagnoses,          & ! [IN]
                                       istep                   ) ! [IN]
          same_mptype_ = .false.
          qnum_flag    = .false.
          uvmet        = .true.
       end select

       if ( rh2qv ) then
          if ( .not. temp2pt ) then
             LOG_ERROR("ParentAtmosInput",*) 'When RH is read, TEMP is necessary'
             call PRC_abort()
          end if

          if ( sfc_diagnoses ) then
             k0 = 2
          else
             k0 = 3
          end if
          !$omp parallel do collapse(2) &
          !$omp private(p_sat,qm)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = k0, KA_org
             if ( temp_org(k,i,j) > UNDEF .and. rh_org(k,i,j) > UNDEF .and. pres_org(k,i,j) > UNDEF ) then
                call psat( temp_org(k,i,j), p_sat )
                qm = EPSvap * rh_org(k,i,j) * 0.01_RP * p_sat &
                     / ( pres_org(k,i,j) - rh_org(k,i,j) * 0.01_RP * p_sat )
                qv_org(k,i,j) = qm / ( 1.0_RP + qm ) ! specific humidity
             else
                qv_org(k,i,j) = UNDEF
             end if
          end do
          end do
          end do
       end if


       select case( upper_qv_type )
       case("COPY")
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 4, KA_org
             if ( qv_org(k,i,j) == UNDEF ) qv_org(k,i,j) = qv_org(k-1,i,j)
          enddo
          enddo
          enddo
       case("ZERO")
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 4, KA_org
             if ( qv_org(k,i,j) == UNDEF ) qv_org(k,i,j) = 0.0_RP
          enddo
          enddo
          enddo
       case default
          LOG_ERROR("ParentAtmosInput",*) 'upper_qv_type in PARAM_MKINIT_REAL is invalid! ', trim(upper_qv_type)
          call PRC_abort
       end select

       if ( mixing_ratio ) then
          !$omp parallel do collapse(2) &
          !$omp private(qtot)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org
             qtot = 0.0_RP
             if ( qv_org(k,i,j) > UNDEF ) then
                qtot = qtot + qv_org(k,i,j)
             end if
             do iq = 1, N_HYD
                if ( qhyd_org(k,i,j,iq) > UNDEF ) then
                   qtot = qtot + qhyd_org(k,i,j,iq)
                end if
             end do
             if( qv_org(k,i,j) > UNDEF ) then
                qv_org(k,i,j) = qv_org(k,i,j) / ( 1.0_RP + qtot )
             else
                qv_org(k,i,j) = 0.0_RP
             end if
             do iq = 1, N_HYD
                if ( qhyd_org(k,i,j,iq) > UNDEF ) then
                   qhyd_org(k,i,j,iq) = qhyd_org(k,i,j,iq) / ( 1.0_RP + qtot )
                else
                   qhyd_org(k,i,j,iq) = 0.0_RP
                end if
             end do
          end do
          end do
          end do
       end if


       if ( .not. same_mptype_ ) then
          if ( qnum_flag ) then
             call ATMOS_PHY_MP_driver_qhyd2qtrc( KA_org, 2, KA_org, IA_org, 1, IA_org, JA_org, 1, JA_org, &
                                                 QV_org(:,:,:), QHYD_org(:,:,:,:), & ! [IN]
                                                 QTRC_org(:,:,:,QS_MP:QE_MP),      & ! [OUT]
                                                 QNUM=QNUM_org(:,:,:,:)            ) ! [IN]
          else
             call ATMOS_PHY_MP_driver_qhyd2qtrc( KA_org, 2, KA_org, IA_org, 1, IA_org, JA_org, 1, JA_org, &
                                                 QV_org(:,:,:), QHYD_org(:,:,:,:), & ! [IN]
                                                 QTRC_org(:,:,:,QS_MP:QE_MP)       ) ! [OUT]
          end if
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             QTRC_org(1,i,j,QS_MP:QE_MP) = UNDEF
             do k = 2, KA_org
                if ( QV_org(k,i,j) == UNDEF ) QTRC_org(k,i,j,QS_MP:QE_MP) = UNDEF
             end do
          end do
          end do
       end if

       if ( pt_dry .and. .not. temp2pt ) then
          if ( nopres ) then
             LOG_ERROR('ParentAtmosInput',*) 'PRES is required'
             call PRC_abort
          end if
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 3, KA_org
             temp_org(k,i,j) = pt_org(k,i,j) * ( pres_org(k,i,j) / PRE00 )**(Rdry/CPdry)
          end do
          end do
          end do
          temp2pt = .true.
       end if

       if ( temp2pt .or. nopres .or. nodens ) then
          call THERMODYN_specific_heat( &
               KA_org, 3, KA_org, IA_org, 1, IA_org, JA_org, 1, JA_org, QA, &
               qtrc_org(:,:,:,:), &
               TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), &
               qdry(:,:,:), Rtot(:,:,:), CVtot(:,:,:), CPtot(:,:,:) )
       end if

       if ( temp2pt ) then
          if ( nopres ) then
             LOG_ERROR('ParentAtmosInput',*) 'If TEMP is read, PRES is required'
             call PRC_abort
          end if
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 3, KA_org
             if ( TEMP_org(k,i,j) == UNDEF ) then
                PT_org(k,i,j) = UNDEF
             else
                call THERMODYN_temp_pres2pott( &
                     TEMP_org(k,i,j), PRES_org(k,i,j), & ! [IN]
                     CPtot(k,i,j), Rtot(k,i,j),        & ! [IN]
                     PT_org(k,i,j)                     ) ! [OUT]
             end if
          enddo
          enddo
          enddo
       endif

       if ( nopres ) then
          if ( nodens ) then
             LOG_ERROR('ParentAtmosInput',*) 'If PRES does not exist, DENS is required'
             call PRC_abort
          end if
          !$omp parallel do collapse(2) &
          !$omp private(rhot_tmp)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 3, KA_org
             if ( PT_org(k,i,j) == UNDEF ) then
                PRES_org(k,i,j) = UNDEF
             else
                rhot_tmp = PT_org(k,i,j) * DENS_org(k,i,j)
                call THERMODYN_rhot2pres( &
                     rhot_tmp, CVtot(k,i,j), CPtot(k,i,j), Rtot(k,i,j), & ! [IN]
                     PRES_org(k,i,j)                                    ) ! [OUT]
             end if
          enddo
          enddo
          enddo
       end if

       if ( nodens .and. use_file_density ) then
          LOG_ERROR('ParentAtmosInput',*) 'DENS is required when USE_FILE_DENSITY is true'
          call PRC_abort
       end if


       lm_layer(:,:) = 3
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 3, KA_org
          ! search the lowermost layer excluding UNDEF
          if( pres_org(k,i,j) == UNDEF ) then
             lm_layer(i,j) = k + 1
          else
             exit
          end if
       end do
       end do
       end do
       if ( sfc_diagnoses ) then

          if ( .not. under_sfc ) then
             !$omp parallel do collapse(2)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = 3, KA_org
                if ( .not. ( &
                     .not. (cz_org(2,i,j) > UNDEF .and. cz_org(k,i,j) > cz_org(2,i,j) ) .or. &
                     .not. (pres_org(2,i,j) > UNDEF .and. pres_org(k,i,j) > pres_org(2,i,j) ) &
                     ) ) then
                   lm_layer(i,j) = k
                   exit
                end if
             end do
             end do
             end do
          end if

          ! groud surface

          !$omp parallel do &
          !$omp private(k,dz)
          do j = 1, JA_org
          do i = 1, IA_org
             k = lm_layer(i,j)

             if ( cz_org(2,i,j) > UNDEF .and. & ! topo exists
                  ( .not. under_sfc .or. cz_org(2,i,j) < cz_org(k,i,j) ) ) then ! surface is lower than the lowest data
                dz = cz_org(k,i,j) - cz_org(2,i,j)
                if ( qv_org(2,i,j) > UNDEF ) qv_org(2,i,j) = qtrc_org(k,i,j,QS_MP)
                if ( temp_org(2,i,j) > UNDEF .and. pres_org(2,i,j) > UNDEF ) then
                   Rtot(2,i,j) = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                   dens_org(2,i,j) = pres_org(2,i,j) / ( Rtot(2,i,j) * temp_org(2,i,j) )
                else if ( pres_org(2,i,j) > UNDEF ) then
                   dens_org(2,i,j) = - ( pres_org(k,i,j) - pres_org(2,i,j) ) * 2.0_RP / ( GRAV * dz ) &
                        - dens_org(k,i,j)
                   Rtot(2,i,j) = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                   temp_org(2,i,j) = pres_org(2,i,j) / ( Rtot(2,i,j) * dens_org(2,i,j) )
                else if ( temp_org(2,i,j) > UNDEF ) then
                   Rtot(2,i,j) = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                   dens_org(2,i,j) = ( pres_org(k,i,j) + GRAV * dens_org(k,i,j) * dz * 0.5_RP ) &
                        / ( Rtot(2,i,j) * temp_org(2,i,j) - GRAV * dz * 0.5_RP )
                   pres_org(2,i,j) = dens_org(2,i,j) * Rtot(2,i,j) * temp_org(2,i,j)
                else
                   temp_org(2,i,j) = temp_org(k,i,j) + LAPS * dz
                   Rtot(2,i,j) = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                   dens_org(2,i,j) = ( pres_org(k,i,j) + GRAV * dens_org(k,i,j) * dz * 0.5_RP ) &
                        / ( Rtot(2,i,j) * temp_org(2,i,j) - GRAV * dz * 0.5_RP )
                   pres_org(2,i,j) = dens_org(2,i,j) * Rtot(2,i,j) * temp_org(2,i,j)
                end if
                CPtot(2,i,j) = CPdry + ( CPvap - CPdry ) * qv_org(2,i,j)
                pt_org(2,i,j) = temp_org(2,i,j) * ( PRE00 / pres_org(2,i,j) )**(Rtot(2,i,j)/CPtot(2,i,j))

             else ! no topo

                ! ignore surface variables
                cz_org  (2,i,j)   = cz_org  (k,i,j)
                w_org   (2,i,j)   = w_org   (k,i,j)
                u_org   (2,i,j)   = u_org   (k,i,j)
                v_org   (2,i,j)   = v_org   (k,i,j)
                pres_org(2,i,j)   = pres_org(k,i,j)
                pt_org  (2,i,j)   = pt_org  (k,i,j)
                dens_org(2,i,j)   = dens_org(k,i,j)
                qtrc_org(2,i,j,:) = qtrc_org(k,i,j,:)

             end if

          end do
          end do


          ! sea level

          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             if ( pres_org(1,i,j) > UNDEF ) then
                temp_org(1,i,j) = temp_org(2,i,j) + LAPS * cz_org(2,i,j)
                Rtot(1,i,j) = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                CPtot(1,i,j) = CPdry + ( CPvap - CPdry ) * qv_org(2,i,j)
                dens_org(1,i,j) = pres_org(1,i,j) / ( Rtot(1,i,j) * temp_org(1,i,j) )
                pt_org(1,i,j) = temp_org(1,i,j) * ( PRE00 / pres_org(1,i,j) )**(Rtot(1,i,j)/CPtot(1,i,j))
                cz_org(1,i,j) = 0.0_RP
                w_org (1,i,j) = w_org(2,i,j)
                u_org (1,i,j) = u_org(2,i,j)
                v_org (1,i,j) = v_org(2,i,j)
                qtrc_org(1,i,j,:) = qtrc_org(2,i,j,:)
             else
                cz_org  (1,i,j)   = UNDEF
                w_org   (1,i,j)   = UNDEF
                u_org   (1,i,j)   = UNDEF
                v_org   (1,i,j)   = UNDEF
                pres_org(1,i,j)   = UNDEF
                pt_org  (1,i,j)   = UNDEF
                dens_org(1,i,j)   = UNDEF
                qtrc_org(1,i,j,:) = UNDEF
             end if
          end do
          end do

       else

          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             cz_org(1:2,i,j)     = UNDEF
             w_org(1:2,i,j)      = UNDEF
             u_org(1:2,i,j)      = UNDEF
             v_org(1:2,i,j)      = UNDEF
             pres_org(1:2,i,j)   = UNDEF
             temp_org(1:2,i,j)   = UNDEF
             dens_org(1:2,i,j)   = UNDEF
             qtrc_org(1:2,i,j,:) = UNDEF
          end do
          end do

       end if ! sfc_diagnoses

    endif ! read by this process?

    call PROF_rapend  ('___AtmosInput',3)

    call PROF_rapstart('___AtmosBcast',3)

    if ( serial_atmos ) then
       if ( first_atmos .or. update_coord ) then
          call COMM_bcast( KA_org, IA_org, JA_org, CZ_org )
       end if
       call COMM_bcast( KA_org, IA_org, JA_org,     PRES_org )
       call COMM_bcast( KA_org, IA_org, JA_org,     W_org    )
       call COMM_bcast( KA_org, IA_org, JA_org,     U_org    )
       call COMM_bcast( KA_org, IA_org, JA_org,     V_org    )
       call COMM_bcast( KA_org, IA_org, JA_org,     PT_org   )
       call COMM_bcast( KA_org, IA_org, JA_org,     DENS_org )
       call COMM_bcast( KA_org, IA_org, JA_org, QA, QTRC_org )
    endif

    call PROF_rapend  ('___AtmosBcast',3)

    !$omp parallel do collapse(3)
    do iq = 1, QA
    do j = 1, JA_org
    do i = 1, IA_org
    do k = 1, KA_org
       if ( QTRC_org(k,i,j,iq) .ne. UNDEF ) then
          QTRC_org(k,i,j,iq) = max( QTRC_org(k,i,j,iq), 0.0_RP )
       end if
    enddo
    enddo
    enddo
    enddo

    ! interpolation
    call PROF_rapstart('___AtmosInterp',3)

    if ( first_atmos .or. update_coord ) then

       k = KA_org
       call INTERP_domain_compatibility( LON_org(:,:),    & ! [IN]
                                         LAT_org(:,:),    & ! [IN]
                                         CZ_org (k,:,:),  & ! [IN]
                                         LON    (:,:),    & ! [IN]
                                         LAT    (:,:),    & ! [IN]
                                         CZ     (KE,:,:), & ! [IN]
                                         FZ     (KE,:,:), & ! [IN]
                                         skip_z = skip_vcheck ) ! [IN]

       select case( itp_type_a )
       case ( i_intrp_linear )

          if ( IA_org == 1 .or. JA_org == 1 ) then
             LOG_ERROR("ParentAtmosInput",*) 'LINER interpolation requires nx, ny > 1'
             LOG_ERROR_CONT(*)               'Use "DIST-WEIGHT" as INTRP_TYPE of PARAM_MKINIT_REAL_ATMOS'
             call PRC_abort
          end if

          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             LAT_org(i,j) = sign( min( abs(LAT_org(i,j)), PI * 0.499999_RP ), LAT_org(i,j) )
          end do
          end do

          call MAPPROJECTION_lonlat2xy( IA_org, 1, IA_org, & ! [IN]
                                        JA_org, 1, JA_org, & ! [IN]
                                        LON_org(:,:),      & ! [IN]
                                        LAT_org(:,:),      & ! [IN]
                                        X_org  (:,:),      & ! [OUT]
                                        Y_org  (:,:)       ) ! [OUT]

          zonal = ( maxval(LON_org) - minval(LON_org) ) > 2.0_RP * PI * 0.9_RP
          pole = ( maxval(LAT_org) > PI * 0.5_RP * 0.9_RP ) .or. ( minval(LAT_org) < - PI * 0.5_RP * 0.9_RP )
          call INTERP_factor3d( KA_org, 1, KA_org,       & ! [IN]
                                IA_org, JA_org,          & ! [IN]
                                KA, KS, KE,              & ! [IN]
                                IA, JA,                  & ! [IN]
                                X_org(:,:), Y_org(:,:),  & ! [IN]
                                CZ_org (:,:,:),          & ! [IN]
                                CX(:), CY(:),            & ! [IN]
                                CZ     (:,:,:),          & ! [IN]
                                igrd   (    :,:,:),      & ! [OUT]
                                jgrd   (    :,:,:),      & ! [OUT]
                                hfact  (    :,:,:),      & ! [OUT]
                                kgrd   (:,:,:,:,:),      & ! [OUT]
                                vfact  (:,  :,:,:),      & ! [OUT]
                                flag_extrap = .false.,   & ! [IN]
                                zonal = zonal,           & ! [IN]
                                pole  = pole             ) ! [IN]

       case ( I_intrp_dstwgt )

          call INTERP_factor3d( itp_nh_a,             & ! [IN]
                                KA_org, 1, KA_org,    & ! [IN]
                                IA_org, JA_org,       & ! [IN]
                                KA, KS, KE,           & ! [IN]
                                IA, JA,               & ! [IN]
                                LON_org(:,:),         & ! [IN]
                                LAT_org(:,:),         & ! [IN]
                                CZ_org (:,:,:),       & ! [IN]
                                LON    (:,:),         & ! [IN]
                                LAT    (:,:),         & ! [IN]
                                CZ     (:,:,:),       & ! [IN]
                                igrd   (    :,:,:),   & ! [OUT]
                                jgrd   (    :,:,:),   & ! [OUT]
                                hfact  (    :,:,:),   & ! [OUT]
                                kgrd   (:,:,:,:,:),   & ! [OUT]
                                vfact  (:,  :,:,:),   & ! [OUT]
                                flag_extrap = .false. ) ! [IN]

       end select

    endif

    call INTERP_interp3d( itp_nh_a,                 &
                          KA_org, 1, KA_org,        &
                          IA_org, JA_org,           &
                          KA, KS, KE,               &
                          IA, JA,                   &
                          igrd(:,:,:), jgrd(:,:,:), & ! [IN]
                          hfact(:,:,:),             & ! [IN]
                          kgrd(:,:,:,:,:),          & ! [IN]
                          vfact(:,:,:,:),           & ! [IN]
                          CZ_org(:,:,:), CZ(:,:,:), & ! [IN]
                          W_org(:,:,:),             & ! [IN]
                          W     (:,:,:),            & ! [OUT]
                          spline = .false.,         & ! [IN]
                          threshold_undef = 1.0_RP, & ! [IN]
                          wsum = wsum(:,:,:),       & ! [OUT]
                          val2 = work(:,:,:)        ) ! [OUT]
    !$omp parallel do collapse(2) &
    !$omp private(kref)
    do j = 1, JA
    do i = 1, IA
!CDIR NOVECTOR
       do k = KS, KA
          if ( W(k,i,j) .ne. UNDEF ) then
             kref = k
             exit
          end if
       end do
       do k = kref-1, KS, -1
          W(k,i,j) = W(k+1,i,j) * log( ( CZ(k,i,j) - topo(i,j) ) / Z0M(i,j) ) / log( ( CZ(k+1,i,j) - topo(i,j) ) / Z0M(i,j) ) * ( 1.0_RP - wsum(k,i,j) ) &
                   + work(k,i,j) * wsum(k,i,j)
       end do
       do k = kref+1, KE
          if ( W(k,i,j) == UNDEF ) W(k,i,j) = W(k-1,i,j)
       end do
    end do
    end do
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                              W(:,:,:), FILTER_ORDER, FILTER_NITER )
       call COMM_vars8( W(:,:,:), 1 )
       call COMM_wait ( W(:,:,:), 1, .false. )
    end if

    if ( .not. uvmet ) then
       ! rotation from map-projected field to latlon field
       !$omp parallel do collapse(2) &
       !$omp private(u_on_map,v_on_map)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org
          if ( U_org(k,i,j) > UNDEF .and. V_org(k,i,j) > UNDEF ) then
             u_on_map =  U_org(k,i,j) * ROTC_org(i,j,1) - V_org(k,i,j) * ROTC_org(i,j,2)
             v_on_map =  U_org(k,i,j) * ROTC_org(i,j,2) + V_org(k,i,j) * ROTC_org(i,j,1)

             U_org(k,i,j) = u_on_map
             V_org(k,i,j) = v_on_map
          end if
       enddo
       enddo
       enddo
    end if


    call INTERP_interp3d( itp_nh_a,                 &
                          KA_org, 1, KA_org,        &
                          IA_org, JA_org,           &
                          KA, KS, KE,               &
                          IA, JA,                   &
                          igrd(:,:,:), jgrd(:,:,:), & ! [IN]
                          hfact(:,:,:),             & ! [IN]
                          kgrd(:,:,:,:,:),          & ! [IN]
                          vfact(:,:,:,:),           & ! [IN]
                          CZ_org(:,:,:), CZ(:,:,:), & ! [IN]
                          U_org(:,:,:),             & ! [IN]
                          U(:,:,:),                 & ! [OUT]
                          spline = .false.,         & ! [IN]
                          threshold_undef = 1.0_RP, & ! [IN]
                          wsum = wsum(:,:,:),       & ! [OUT]
                          val2 = work(:,:,:)        ) ! [OUT]
    !$omp parallel do collapse(2) &
    !$omp private(kref)
    do j = 1, JA
    do i = 1, IA
       do k = KS, KA
          if ( U(k,i,j) .ne. UNDEF ) then
             kref = k
             exit
          end if
       end do
       do k = kref-1, KS, -1
          U(k,i,j) = U(k+1,i,j) * log( ( CZ(k,i,j) - topo(i,j) ) / Z0M(i,j) ) / log( ( CZ(k+1,i,j) - topo(i,j) ) / Z0M(i,j) ) * ( 1.0_RP - wsum(k,i,j) ) &
                   + work(k,i,j) * wsum(k,i,j)
       end do
       do k = kref+1, KE
          if ( U(k,i,j) == UNDEF ) U(k,i,j) = U(k-1,i,j)
       end do
    end do
    end do
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                              U(:,:,:), FILTER_ORDER, FILTER_NITER )
       call COMM_vars8( U(:,:,:), 1 )
       call COMM_wait ( U(:,:,:), 1, .false. )
    end if

    call INTERP_interp3d( itp_nh_a,                 &
                          KA_org, 1, KA_org,        &
                          IA_org, JA_org,           &
                          KA, KS, KE,               &
                          IA, JA,                   &
                          igrd(:,:,:), jgrd(:,:,:), & ! [IN]
                          hfact(:,:,:),             & ! [IN]
                          kgrd(:,:,:,:,:),          & ! [IN]
                          vfact(:,:,:,:),           & ! [IN]
                          CZ_org(:,:,:), CZ(:,:,:), & ! [IN]
                          V_org(:,:,:),             & ! [IN]
                          V(:,:,:),                 & ! [OUT]
                          spline = .false.,         & ! [IN]
                          threshold_undef = 1.0_RP, & ! [IN]
                          wsum = wsum(:,:,:),       & ! [OUT]
                          val2 = work(:,:,:)        ) ! [OUT]
    !$omp parallel do collapse(2) &
    !$omp private(kref)
    do j = 1, JA
    do i = 1, IA
       do k = KS, KA
          if ( V(k,i,j) .ne. UNDEF ) then
             kref = k
             exit
          end if
       end do
       do k = kref-1, KS, -1
          V(k,i,j) = V(k+1,i,j) * log( ( CZ(k,i,j) - topo(i,j) ) / Z0M(i,j) ) / log( ( CZ(k+1,i,j) - topo(i,j) ) / Z0M(i,j) ) * ( 1.0_RP - wsum(k,i,j) ) &
                   + work(k,i,j) * wsum(k,i,j)
       end do
       do k = kref+1, KE
          if ( V(k,i,j) == UNDEF ) V(k,i,j) = V(k-1,i,j)
       end do
    end do
    end do
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                              V(:,:,:), FILTER_ORDER, FILTER_NITER )
       call COMM_vars8( V(:,:,:), 1 )
       call COMM_wait ( V(:,:,:), 1, .false. )
    end if

    ! rotation from latlon field to map-projected field
    !$omp parallel do collapse(2) &
    !$omp private(u_on_map,v_on_map)
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

    ! from scalar point to staggered point
    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE-1
       VELZ(k,i,j) = 0.5_RP * ( W(k+1,i,j) + W(k,i,j) )
    enddo
    enddo
    enddo

    !$omp parallel do
    do j = 1, JA
    do i = 1, IA-1
    do k = KS, KE
       VELX(k,i,j) = 0.5_RP * ( U(k,i+1,j) + U(k,i,j) )
    enddo
    enddo
    enddo

    i = IA
    !$omp parallel do
    do j = 1, JA
    do k = KS, KE
       VELX(k,i,j) = U(k,i,j)
    enddo
    enddo

    !$omp parallel do
    do j = 1, JA-1
    do i = 1, IA
    do k = KS, KE
       VELY(k,i,j) = 0.5_RP * ( V(k,i,j+1) + V(k,i,j) )
    enddo
    enddo
    enddo

    j = JA
    !$omp parallel do
    do i = 1, IA
    do k = KS, KE
       VELY(k,i,j) = V(k,i,j)
    enddo
    enddo

    !$omp parallel do
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

    call INTERP_interp3d( itp_nh_a,                 &
                          KA_org, 1, KA_org,        &
                          IA_org, JA_org,           &
                          KA, KS, KE,               &
                          IA, JA,                   &
                          igrd(:,:,:), jgrd(:,:,:), & ! [IN]
                          hfact(:,:,:),             & ! [IN]
                          kgrd(:,:,:,:,:),          & ! [IN]
                          vfact(:,:,:,:),           & ! [IN]
                          CZ_org(:,:,:), CZ(:,:,:), & ! [IN]
                          PT_org(:,:,:),            & ! [IN]
                          POTT(:,:,:),              & ! [OUT]
                          spline = .false.,         & ! [IN]
                          threshold_undef = 1.0_RP, & ! [IN]
                          wsum = wsum(:,:,:),       & ! [OUT]
                          val2 = work(:,:,:)        ) ! [OUT]
    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
       do k = KS+1, KE
          if ( POTT(k,i,j) == UNDEF .and. POTT(k-1,i,j) .ne. UNDEF ) POTT(k,i,j) = POTT(k-1,i,j)
       end do
       do k = KE-1, KS, -1
          if ( POTT(k,i,j) == UNDEF ) then
             POTT(k,i,j) = POTT(k+1,i,j) * ( 1.0_RP - wsum(k,i,j) ) &
                         + work(k,i,j) * wsum(k,i,j)
          end if
       end do
       POTT(   1:KS-1,i,j) = UNDEF
       POTT(KE+1:KA  ,i,j) = UNDEF
    enddo
    enddo
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                              POTT(:,:,:), FILTER_ORDER, FILTER_NITER )
       call COMM_vars8( POTT(:,:,:), 1 )
       call COMM_wait ( POTT(:,:,:), 1, .false. )
    end if

    do iq = 1, QA

       if ( ( iq < QS_MP .or. iq > QE_MP ) .and. ( .not. qtrc_flag(iq) ) ) then
          !$omp parallel do collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = 1, KA
             QTRC(k,i,j,iq) = UNDEF
          end do
          end do
          end do
          cycle
       end if

       call INTERP_interp3d( itp_nh_a,                 &
                             KA_org, 1, KA_org,        &
                             IA_org, JA_org,           &
                             KA, KS, KE,               &
                             IA, JA,                   &
                             igrd(:,:,:), jgrd(:,:,:), & ! [IN]
                             hfact(:,:,:),             & ! [IN]
                             kgrd(:,:,:,:,:),          & ! [IN]
                             vfact(:,:,:,:),           & ! [IN]
                             CZ_org(:,:,:), CZ(:,:,:), & ! [IN]
                             QTRC_org(:,:,:,iq),       & ! [IN]
                             QTRC(:,:,:,iq),           & ! [OUT]
                             spline = .false.,         & ! [IN]
                             threshold_undef = 1.0_RP, & ! [IN]
                             wsum = wsum(:,:,:),       & ! [OUT]
                             val2 = work(:,:,:)        ) ! [OUT]
       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
          do k = KS+1, KE
             if ( QTRC(k,i,j,iq) == UNDEF .and. QTRC(k-1,i,j,iq) > UNDEF ) QTRC(k,i,j,iq) = QTRC(k-1,i,j,iq)
          end do
          do k = KE-1, KS, -1
             if ( QTRC(k,i,j,iq) == UNDEF .and. QTRC(k+1,i,j,iq) > UNDEF ) then
                QTRC(k,i,j,iq) = QTRC(k+1,i,j,iq) * ( 1.0_RP - wsum(k,i,j) ) &
                               + work(k,i,j) * wsum(k,i,j)
             end if
          end do
          do k = KS, KE
             QTRC(k,i,j,iq) = max( QTRC(k,i,j,iq), 0.0_RP )
          end do
          QTRC(   1:KS-1,i,j,iq) = 0.0_RP
          QTRC(KE+1:KA  ,i,j,iq) = 0.0_RP
       enddo
       enddo
       if ( FILTER_NITER > 0 ) then
          !$omp parallel do collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = 1, KA
             one(k,i,j) = 1.0_RP
          end do
          end do
          end do
          call FILTER_hyperdiff( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                 QTRC(:,:,:,iq), FILTER_ORDER, FILTER_NITER, &
                                 limiter_sign = one(:,:,:) )
          call COMM_vars8( QTRC(:,:,:,iq), 1 )
          call COMM_wait ( QTRC(:,:,:,iq), 1, .false. )
       end if
    enddo

    call INTERP_interp3d( itp_nh_a,            &
                          KA_org, 1, KA_org,   &
                          IA_org, JA_org,      &
                          KA, KS, KE,          &
                          IA, JA,              &
                          igrd    (    :,:,:), & ! [IN]
                          jgrd    (    :,:,:), & ! [IN]
                          hfact   (    :,:,:), & ! [IN]
                          kgrd    (:,:,:,:,:), & ! [IN]
                          vfact   (:,  :,:,:), & ! [IN]
                          CZ_org  (:,:,:),     & ! [IN]
                          CZ      (:,:,:),     & ! [IN]
                          PRES_org(:,:,:),     & ! [IN]
                          PRES    (:,:,:),     & ! [OUT]
                          logwgt = .true.      ) ! [IN, optional]

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       QC(k,i,j) = 0.0_RP
    end do
    end do
    end do
    if ( ATMOS_HYDROMETEOR_dry ) then
       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          QV(k,i,j) = 0.0_RP
       end do
       end do
       end do
    else
       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          QV(k,i,j) = QTRC(k,i,j,I_QV)
          do iq = QLS, QLE
             QC(k,i,j) = QC(k,i,j) + QTRC(k,i,j,iq)
          enddo
       end do
       end do
       end do
    end if


    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       PRES2(k,i,j) = PRES(k,i,j)
    end do
    end do
    end do

    if ( use_file_density ) then
       call INTERP_interp3d( itp_nh_a,                 &
                             KA_org, 1, KA_org,        &
                             IA_org, JA_org,           &
                             KA, KS, KE,               &
                             IA, JA,                   &
                             igrd(:,:,:), jgrd(:,:,:), & ! [IN]
                             hfact(:,:,:),             & ! [IN]
                             kgrd(:,:,:,:,:),          & ! [IN]
                             vfact(:,:,:,:),           & ! [IN]
                             CZ_org(:,:,:), CZ(:,:,:), & ! [IN]
                             DENS_org(:,:,:),          & ! [IN]
                             DENS(:,:,:),              & ! [OUT]
                             threshold_undef = 1.0_RP, & ! [IN]
                             wsum = wsum(:,:,:),       & ! [OUT]
                             val2 = work(:,:,:)        ) ! [OUT]
       call HYDROSTATIC_buildrho_real( KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                                       POTT(:,:,:), QV(:,:,:), QC(:,:,:), & ! [IN]
                                       CZ(:,:,:),                         & ! [IN]
                                       PRES2(:,:,:),                      & ! [INOUT]
                                       DENS2(:,:,:), TEMP(:,:,:)          ) ! [OUT]
       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          if ( DENS(k,i,j) == UNDEF ) then
             DENS(k,i,j) = DENS2(k,i,j) * ( 1.0_RP - wsum(k,i,j) ) &
                         + work(k,i,j) * wsum(k,i,j)
          end if
       end do
       end do
       end do
    else
       ! make density & pressure profile in moist condition
       call HYDROSTATIC_buildrho_real( KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                                       POTT(:,:,:), QV(:,:,:), QC(:,:,:), & ! [IN]
                                       CZ(:,:,:),                         & ! [IN]
                                       PRES2(:,:,:),                      & ! [INOUT]
                                       DENS(:,:,:), TEMP(:,:,:)           ) ! [OUT]
    endif

    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                              DENS(:,:,:), FILTER_ORDER, FILTER_NITER, &
                              limiter_sign = one(:,:,:) )
       call COMM_vars8( DENS(:,:,:), 1 )
       call COMM_wait ( DENS(:,:,:), 1, .false. )
    end if

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       if ( PRES(k,i,j) == UNDEF ) PRES(k,i,j) = PRES2(k,i,j)
    end do
    end do
    end do


    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       DENS(   1:KS-1,i,j) = 0.0_RP
       DENS(KE+1:KA  ,i,j) = 0.0_RP
    enddo
    enddo

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE-1
       MOMZ(k,i,j) = VELZ(k,i,j) * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    !$omp parallel do
    do j = 1, JA
    do i = 1, IA-1
    do k = KS, KE
       MOMX(k,i,j) = VELX(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    i = IA
    !$omp parallel do
    do j = 1, JA
    do k = KS, KE
       MOMX(k,i,j) = VELX(k,i,j) * DENS(k,i,j)
    enddo
    enddo

    !$omp parallel do
    do j = 1, JA-1
    do i = 1, IA
    do k = KS, KE
       MOMY(k,i,j) = VELY(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    j = JA
    !$omp parallel do
    do i = 1, IA
    do k = KS, KE
       MOMY(k,i,j) = VELY(k,i,j) * DENS(k,i,j)
    enddo
    enddo

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RHOT(k,i,j) = POTT(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do
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

    first_atmos = .false.

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
       qtrc_flag, &
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
    logical,          intent(in)  :: qtrc_flag(QA)
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
         'PT',   'Reference PT',      'K',     'ZXYT',  datatype, & ! [IN]
         vid(5),                                                  & ! [OUT]
         timeintv=timeintv                                        ) ! [IN]

    do iq = QS_MP, QE_MP
       call FILE_CARTESC_def_var( fid,                               & ! [IN]
            TRACER_NAME(iq), 'Reference '//TRACER_NAME(iq), 'kg/kg', & ! [IN]
            'ZXYT', datatype,                                        & ! [IN]
            vid(5+iq),                                               & ! [OUT]
            timeintv = timeintv                                      ) ! [IN]
    enddo

    do iq = 1, QA
       if ( iq >= QS_MP .and. iq <= QE_MP ) cycle
       if ( .not. qtrc_flag(iq) ) cycle
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
       DENS,      &
       VELZ,      &
       VELX,      &
       VELY,      &
       POTT,      &
       QTRC,      &
       qtrc_flag, &
       fid, vid,  &
       timeintv,  &
       istep      )
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: VELZ(KA,IA,JA)
    real(RP), intent(in) :: VELX(KA,IA,JA)
    real(RP), intent(in) :: VELY(KA,IA,JA)
    real(RP), intent(in) :: POTT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    logical,  intent(in) :: qtrc_flag(QA)
    integer,  intent(in) :: fid
    integer,  intent(in) :: vid(5+QA)
    real(DP), intent(in) :: timeintv
    integer,  intent(in) :: istep

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
    call FILE_CARTESC_write_var( fid, vid(5), work(:,:,:,:), 'PT',   'ZXYT', timeintv, timeofs=timeofs )

    do iq = QS_MP, QE_MP
       call FILE_CARTESC_write_var( fid, vid(5+iq),QTRC(:,:,:,iq:iq), TRACER_NAME(iq), &
                              'ZXYT', timeintv, timeofs=timeofs                  )
    enddo

    do iq = 1, QA
       if ( iq >= QS_MP .and. iq <= QS_MP ) cycle
       if ( .not. qtrc_flag(iq) ) cycle
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
       basename_org,        &
       basename_land,       &
       basename_ocean,      &
       filetype_land,       &
       filetype_ocean,      &
       use_file_landwater,  &
       intrp_iter_max,      &
       intrp_land_temp,     &
       intrp_land_water,    &
       intrp_land_sfc_temp, &
       intrp_ocean_temp,    &
       intrp_ocean_sfc_temp )
    use scale_comm_cartesC, only: &
       COMM_world, &
       COMM_datatype
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LON, &
       ATMOS_GRID_CARTESC_REAL_LAT
    use mod_realinput_netcdf, only: &
         ParentLandSetupNetCDF, &
         ParentOceanSetupNetCDF
    use mod_realinput_grads, only: &
         ParentLandSetupGrADS, &
         ParentOceanSetupGrADS
    implicit none

    integer,          intent(out) :: ldims(3) ! dims for land
    integer,          intent(out) :: odims(2) ! dims for ocean
    integer,          intent(out) :: lmdlid   ! model id for land
    integer,          intent(out) :: omdlid   ! model id for ocean
    integer,          intent(out) :: timelen  ! number of time steps in ocean file

    character(len=*), intent(in)  :: basename_org
    character(len=*), intent(in)  :: basename_land
    character(len=*), intent(in)  :: basename_ocean
    character(len=*), intent(in)  :: filetype_land
    character(len=*), intent(in)  :: filetype_ocean
    logical,          intent(in)  :: use_file_landwater ! use land water data from files
    integer,          intent(in)  :: intrp_iter_max
    character(len=*), intent(in)  :: intrp_land_temp
    character(len=*), intent(in)  :: intrp_land_water
    character(len=*), intent(in)  :: intrp_land_sfc_temp
    character(len=*), intent(in)  :: intrp_ocean_temp
    character(len=*), intent(in)  :: intrp_ocean_sfc_temp

    real(RP), pointer :: lon_all(:,:)
    real(RP), pointer :: lat_all(:,:)

    real(RP) :: LON_min, LON_max
    real(RP) :: LAT_min, LAT_max

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentSurfaceSetup",*) 'Setup'

    ! Land

    if( LKMAX < 4 )then
       LOG_ERROR("ParentSurfaceSetup",*) 'LKMAX less than 4: ', LKMAX
       LOG_ERROR_CONT(*) 'in Real Case, LKMAX should be set more than 4'
       call PRC_abort
    endif

    LOG_INFO("ParentSurfaceSetup",*) 'Horizontal Interpolation Level: ', COMM_CARTESC_NEST_INTERP_LEVEL


    if( serial_land .and. .not. PRC_IsMaster ) then
       do_read_land = .false.
    else
       do_read_land = .true.
    endif

    select case(filetype_land)
    case('NetCDF')

       lmdlid = iNetCDF
       call ParentLandSetupNetCDF( ldims,              & ! (out)
                                   timelen,            & ! (out)
                                   lon_all, lat_all,   & ! (out)
                                   basename_org,       & ! (in)
                                   basename_land,      & ! (in)
                                   use_file_landwater, & ! (in)
                                   serial_land,        & ! (inout)
                                   do_read_land        ) ! (inout)

    case('GrADS')

       lmdlid = iGrADS
       if ( do_read_land ) call ParentLandSetupGrADS( ldims,             & ! (out)
                                                      timelen,           & ! (out)
                                                      lon_all, lat_all,  & ! (out)
                                                      basename_org,      & ! (in)
                                                      basename_land      ) ! (in)

    case default

       LOG_ERROR("ParentSurfaceSetup",*) 'Unsupported FILE TYPE:', trim(filetype_land)
       call PRC_abort

    endselect

    LON_min = minval( ATMOS_GRID_CARTESC_REAL_LON(:,:) )
    LON_max = maxval( ATMOS_GRID_CARTESC_REAL_LON(:,:) )
    LAT_min = minval( ATMOS_GRID_CARTESC_REAL_LAT(:,:) )
    LAT_max = maxval( ATMOS_GRID_CARTESC_REAL_LAT(:,:) )

    if ( serial_land ) then
       call COMM_bcast( 3, ldims(:) )
       call COMM_bcast( timelen )

       if ( .not. do_read_land ) then
          allocate( LON_all(ldims(2), ldims(3)) )
          allocate( LAT_all(ldims(2), ldims(3)) )
       end if
       call COMM_bcast( ldims(2), ldims(3), LON_all )
       call COMM_bcast( ldims(2), ldims(3), LAT_all )

       call MPI_AllReduce( MPI_IN_PLACE, LON_min, 1, COMM_datatype, MPI_MIN, COMM_WORLD, ierr )
       call MPI_AllReduce( MPI_IN_PLACE, LON_max, 1, COMM_datatype, MPI_MAX, COMM_WORLD, ierr )
       call MPI_ALLReduce( MPI_IN_PLACE, LAT_min, 1, COMM_datatype, MPI_MIN, COMM_WORLD, ierr )
       call MPI_ALLReduce( MPI_IN_PLACE, LAT_max, 1, COMM_datatype, MPI_MAX, COMM_WORLD, ierr )
    endif

    call get_IJrange( lIS_org, lIE_org, lJS_org, lJE_org, & ! [OUT]
                      ldims(2), ldims(3),                 & ! [IN]
                      LON_min, LON_max, LAT_min, LAT_max, & ! [IN]
                      LON_all, LAT_all                    ) ! [IN]
    lIS_org = max( lIS_org - INTRP_ITER_MAX, 1 )
    lIE_org = min( lIE_org + INTRP_ITER_MAX, ldims(2) )
    lJS_org = max( lJS_org - INTRP_ITER_MAX, 1 )
    lJE_org = min( lJE_org + INTRP_ITER_MAX, ldims(3) )

    lKS_org = 1
    lKE_org = ldims(1)
    lKA_org = lKE_org - lKS_org + 1
    lIA_org = lIE_org - lIS_org + 1
    lJA_org = lJE_org - lJS_org + 1

    allocate( lLON_org(lIA_org, lJA_org) )
    allocate( lLAT_org(lIA_org, lJA_org) )

    !$omp parallel do
    do j = 1, lJA_org
    do i = 1, lIA_org
       lLON_org(i,j) = LON_all(i-1+lIS_org,j-1+lJS_org)
       lLAT_org(i,j) = LAT_all(i-1+lIS_org,j-1+lJS_org)
    end do
    end do

    deallocate( LON_all )
    deallocate( LAT_all )


    select case( INTRP_LAND_TEMP )
    case( 'off' )
       i_intrp_land_temp = i_intrp_off
    case( 'mask' )
       i_intrp_land_temp = i_intrp_mask
    case( 'fill' )
       i_intrp_land_temp = i_intrp_fill
    case default
       LOG_ERROR("ParentSurfaceSetup",*) 'INTRP_LAND_TEMP is invalid. ', INTRP_LAND_TEMP
       call PRC_abort
    end select
    select case( INTRP_LAND_SFC_TEMP )
    case( 'off' )
       i_intrp_land_sfc_temp = i_intrp_off
    case( 'mask' )
       i_intrp_land_sfc_temp = i_intrp_mask
    case( 'fill' )
       i_intrp_land_sfc_temp = i_intrp_fill
    case default
       LOG_ERROR("ParentSurfaceSetup",*) 'INTRP_LAND_SFC_TEMP is invalid. ', INTRP_LAND_SFC_TEMP
       call PRC_abort
    end select
    select case( INTRP_LAND_WATER )
    case( 'off' )
       i_intrp_land_water = i_intrp_off
    case( 'mask' )
       i_intrp_land_water = i_intrp_mask
    case( 'fill' )
       i_intrp_land_water = i_intrp_fill
    case default
       LOG_ERROR("ParentSurfaceSetup",*) 'INTRP_LAND_WATER is invalid. ', INTRP_LAND_WATER
       call PRC_abort
    end select

    select case( lmdlid )
    case( iNetCDF )
       i_intrp_land_temp      = i_intrp_mask
       i_intrp_land_sfc_temp  = i_intrp_mask
       i_intrp_land_water     = i_intrp_mask
    end select


    ! Ocean

    if( serial_ocean .and. .not. PRC_IsMaster ) then
       do_read_ocean = .false.
    else
       do_read_ocean = .true.
    endif

    select case(filetype_ocean)
    case('NetCDF')

       omdlid = iNetCDF
       call ParentOceanSetupNetCDF( odims, timelen,   & ! (out)
                                    lon_all, lat_all, & ! (out)
                                    basename_org,     & ! (in)
                                    basename_ocean,   & ! (in)
                                    serial_ocean,     & ! (inout)
                                    do_read_ocean     ) ! (inout)

    case('GrADS')

       omdlid = iGrADS
       if ( do_read_ocean ) call ParentOceanSetupGrADS( odims, timelen,   & ! (out)
                                                        lon_all, lat_all, & ! (out)
                                                        basename_org,     & ! (in)
                                                        basename_ocean    ) ! (in)
    case default

       LOG_ERROR("ParentSurfaceSetup",*) 'Unsupported FILE TYPE:', trim(filetype_ocean)
       call PRC_abort

    endselect

    LON_min = minval( ATMOS_GRID_CARTESC_REAL_LON(:,:) )
    LON_max = maxval( ATMOS_GRID_CARTESC_REAL_LON(:,:) )
    LAT_min = minval( ATMOS_GRID_CARTESC_REAL_LAT(:,:) )
    LAT_max = maxval( ATMOS_GRID_CARTESC_REAL_LAT(:,:) )

    if ( serial_ocean ) then
       call COMM_bcast( 2, odims(:) )
       call COMM_bcast( timelen )

       if ( .not. do_read_ocean ) then
          allocate( LON_all(odims(1), odims(2)) )
          allocate( LAT_all(odims(1), odims(2)) )
       end if
       call COMM_bcast( odims(1), odims(2), LON_all )
       call COMM_bcast( odims(1), odims(2), LAT_all )

       call MPI_AllReduce( MPI_IN_PLACE, LON_min, 1, COMM_datatype, MPI_MIN, COMM_WORLD, ierr )
       call MPI_AllReduce( MPI_IN_PLACE, LON_max, 1, COMM_datatype, MPI_MAX, COMM_WORLD, ierr )
       call MPI_ALLReduce( MPI_IN_PLACE, LAT_min, 1, COMM_datatype, MPI_MIN, COMM_WORLD, ierr )
       call MPI_ALLReduce( MPI_IN_PLACE, LAT_max, 1, COMM_datatype, MPI_MAX, COMM_WORLD, ierr )
    endif

    call get_IJrange( oIS_org, oIE_org, oJS_org, oJE_org, & ! [OUT]
                      odims(1), odims(2),                 & ! [IN]
                      LON_min, LON_max, LAT_min, LAT_max, & ! [IN]
                      LON_all, LAT_all                    ) ! [IN]
    oIS_org = max( oIS_org - INTRP_ITER_MAX, 1 )
    oIE_org = min( oIE_org + INTRP_ITER_MAX, odims(1) )
    oJS_org = max( oJS_org - INTRP_ITER_MAX, 1 )
    oJE_org = min( oJE_org + INTRP_ITER_MAX, odims(2) )

    oIA_org = oIE_org - oIS_org + 1
    oJA_org = oJE_org - oJS_org + 1

    allocate( oLON_org(oIA_org, oJA_org) )
    allocate( oLAT_org(oIA_org, oJA_org) )

    !$omp parallel do
    do j = 1, oJA_org
    do i = 1, oIA_org
       oLON_org(i,j) = LON_all(i-1+oIS_org,j-1+oJS_org)
       oLAT_org(i,j) = LAT_all(i-1+oIS_org,j-1+oJS_org)
    end do
    end do

    deallocate( LON_all )
    deallocate( LAT_all )

    select case( INTRP_OCEAN_TEMP )
    case( 'off' )
       i_intrp_ocean_temp = i_intrp_off
    case( 'mask' )
       i_intrp_ocean_temp = i_intrp_mask
    case( 'fill' )
       i_intrp_ocean_temp = i_intrp_fill
    case default
       LOG_ERROR("ParentSurfaceSetup",*) 'INTRP_OCEAN_TEMP is invalid. ', INTRP_OCEAN_TEMP
       call PRC_abort
    end select
    select case( INTRP_OCEAN_SFC_TEMP )
    case( 'off' )
       i_intrp_ocean_sfc_temp = i_intrp_off
    case( 'mask' )
       i_intrp_ocean_sfc_temp = i_intrp_mask
    case( 'fill' )
       i_intrp_ocean_sfc_temp = i_intrp_fill
    case default
       LOG_ERROR("ParentSurfaceSetup",*) 'INTRP_OCEAN_SFC_TEMP is invalid. ', INTRP_OCEAN_SFC_TEMP
       call PRC_abort
    end select

    select case( omdlid )
    case( iNetCDF )
       i_intrp_ocean_temp     = i_intrp_mask
       i_intrp_ocean_sfc_temp = i_intrp_mask
    end select


    allocate( oigrd (       IA, JA, itp_nh_a) )
    allocate( ojgrd (       IA, JA, itp_nh_a) )
    allocate( ohfact(       IA, JA, itp_nh_a) )

    allocate( hfact_ol(oIA_org, oJA_org, itp_nh_ol) )
    allocate( igrd_ol (oIA_org, oJA_org, itp_nh_ol) )
    allocate( jgrd_ol (oIA_org, oJA_org, itp_nh_ol) )

    return
  end subroutine ParentSurfaceSetup

  !-----------------------------------------------------------------------------
  !> Surface Open
  subroutine ParentSurfaceOpen( &
       filetype_land, filetype_ocean, &
       basename_org, basename_land, basename_ocean )
    use mod_realinput_netcdf, only: &
       ParentLandOpenNetCDF, &
       ParentOceanOpenNetCDF
    implicit none
    character(len=*), intent(in) :: filetype_land
    character(len=*), intent(in) :: filetype_ocean
    character(len=*), intent(in) :: basename_org
    character(len=*), intent(in) :: basename_land
    character(len=*), intent(in) :: basename_ocean

    select case ( filetype_land )
    case ( "NetCDF" )
       call ParentLandOpenNetCDF( basename_org, basename_land )
    case ( "GrADS" )
       ! do nothing
    end select

    select case ( filetype_ocean )
    case ( "NetCDF" )
       call ParentOceanOpenNetCDF( basename_org, basename_land )
    case ( "GrADS" )
       ! do nothing
    end select

    return
  end subroutine ParentSurfaceOpen
  !-----------------------------------------------------------------------------
  !> Surface Finalize
  subroutine ParentSurfaceFinalize( &
       filetype_land,               &
       filetype_ocean               )
    use mod_realinput_netcdf, only: &
         ParentLandFinalizeNetCDF, &
         ParentOceanFinalizeNetCDF
    implicit none

    character(len=*), intent(in)  :: filetype_land
    character(len=*), intent(in)  :: filetype_ocean
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentSurfaceFinalize",*) 'Finalize'

    ! Land

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
    case('NetCDF')

       if ( do_read_ocean ) then
          call ParentLandFinalizeNetCDF
       end if

    case('GrADS')

       ! do nothing

    case default

       LOG_ERROR("ParentSurfaceFinalize",*) 'Unsupported FILE TYPE:', trim(filetype_land)
       call PRC_abort

    endselect


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
    case('NetCDF')

       if ( do_read_ocean ) then
          call ParentOceanFinalizeNetCDF
       end if

    case('GrADS')

       ! do nothing

    case default

       LOG_ERROR("ParentSurfaceFinalize",*) 'Unsupported FILE TYPE:', trim(filetype_ocean)
       call PRC_abort

    endselect


    deallocate( llon_org )
    deallocate( llat_org )

    deallocate( olon_org )
    deallocate( olat_org )

    deallocate( oigrd )
    deallocate( ojgrd )
    deallocate( ohfact )

    deallocate( hfact_ol )
    deallocate( igrd_ol  )
    deallocate( jgrd_ol  )

    first_surface = .true.

    return
  end subroutine ParentSurfaceFinalize

  !-----------------------------------------------------------------------------
  !> Boundary Data Write
  subroutine BoundarySurfaceSetup( &
       basename,    &
       title,       &
       timeintv,    &
       multi_ocean, &
       multi_land,  &
       fid,         &
       vid          )
    use scale_file_cartesC, only: &
       FILE_CARTESC_create, &
       FILE_CARTESC_def_var, &
       FILE_CARTESC_enddef
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE
    implicit none
    character(len=*), intent(in) :: basename
    character(len=*), intent(in) :: title
    real(DP),         intent(in) :: timeintv
    logical,          intent(in) :: multi_ocean
    logical,          intent(in) :: multi_land
    integer,          intent(out) :: fid
    integer,          intent(out) :: vid(10)

    character(len=H_SHORT) :: boundary_out_dtype = 'DEFAULT'  !< REAL4 or REAL8
    !---------------------------------------------------------------------------

    call FILE_CARTESC_create( basename, title, boundary_out_dtype, fid, date=nowdate )

    if ( multi_land ) then
       call FILE_CARTESC_def_var( fid,                      & ! [IN]
            'LAND_TEMP', 'Reference Land Temperature', 'K', & ! [IN]
            'LXYT', boundary_out_dtype,                     & ! [IN]
            vid(1),                                         & ! [OUT]
            timeintv=timeintv                               ) ! [IN]
       call FILE_CARTESC_def_var( fid,                        & ! [IN]
            'LAND_WATER', 'Reference Land Moisture', 'm3/m3', & ! [IN]
            'LXYT', boundary_out_dtype,                       & ! [IN]
            vid(2),                                           & ! [OUT]
            timeintv=timeintv                                 ) ! [IN]
       call FILE_CARTESC_def_var( fid,                                  & ! [IN]
            'LAND_SFC_TEMP', 'Reference Land Surface Temperature', 'K', & ! [IN]
            'XYT', boundary_out_dtype,                                  & ! [IN]
            vid(3),                                                     & ! [OUT]
            timeintv=timeintv                                           ) ! [IN]
    end if

    if ( multi_ocean ) then
       call FILE_CARTESC_def_var( fid,                        & ! [IN]
            'OCEAN_TEMP', 'Reference Ocean Temperature', 'K', & ! [IN]
            'OXYT', boundary_out_dtype,                      & ! [IN]
            vid(6),                                           & ! [OUT]
            timeintv=timeintv                                 ) ! [IN]
       call FILE_CARTESC_def_var( fid,                                    & ! [IN]
            'OCEAN_SFC_TEMP', 'Reference Ocean Surface Temperature', 'K', & ! [IN]
            'XYT', boundary_out_dtype,                                   & ! [IN]
            vid(7),                                                       & ! [OUT]
            timeintv=timeintv                                             ) ! [IN]
       call FILE_CARTESC_def_var( fid,                         & ! [IN]
            'OCEAN_SFC_Z0', 'Reference Ocean Surface Z0', 'm', & ! [IN]
            'XYT', boundary_out_dtype,                        & ! [IN]
            vid(10),                                           & ! [OUT]
            timeintv=timeintv                                  ) ! [IN]
    end if

    call FILE_CARTESC_enddef( fid )

    return
  end subroutine BoundarySurfaceSetup

  !-----------------------------------------------------------------------------
  !> Surface Data Read
  subroutine ParentSurfaceInput( &
       tg, strg, lst, albg,               &
       tc_urb, qc_urb, uc_urb, ust, albu, &
       lst_ocean,                         &
       lz_org, topo_org,                  &
       lmask_org, omask_org,              &
       tw, sst, albw, z0w,                &
       basename_org,                      &
       basename_land, basename_ocean,     &
       mdlid_land, mdlid_ocean,           &
       ldims, odims,                      &
       use_file_landwater,                &
       init_landwater_ratio,              &
!       init_landwater_ratio_each,         &
       init_ocean_alb_lw,                 &
       init_ocean_alb_sw,                 &
       init_ocean_z0w,                    &
       intrp_iter_max,                    &
       soilwater_ds2vc_flag,              &
       elevation_correction_land,         &
       elevation_correction_ocean,        &
       oistep, listep,                    &
       multi_land,                        &
       URBAN_do                           )
    use scale_comm_cartesC, only: &
         COMM_bcast, &
         COMM_vars8, &
         COMM_wait
    use scale_const, only: &
         UNDEF => CONST_UNDEF, &
         EPS   => CONST_EPS, &
         LAPS  => CONST_LAPS
    use scale_topography, only: &
         TOPOGRAPHY_Zsfc
    use scale_interp, only: &
         INTERP_factor2d, &
         INTERP_interp2d
    use scale_atmos_grid_cartesC, only: &
         CX  => ATMOS_GRID_CARTESC_CX, &
         CY  => ATMOS_GRID_CARTESC_CY
    use scale_land_grid_cartesC, only: &
         LCZ => LAND_GRID_CARTESC_CZ
    use scale_landuse, only: &
         lsmask_nest => LANDUSE_frac_land, &
         LANDUSE_PFT_nmax
    use mod_realinput_netcdf, only: &
         ParentOCeanInputNetCDF, &
         ParentLandInputNetCDF
    use mod_realinput_grads, only: &
         ParentOceanInputGrADS, &
         ParentLandInputGrADS
    implicit none

    real(RP),         intent(inout) :: tg  (LKA,IA,JA)
    real(RP),         intent(inout) :: strg(LKA,IA,JA)
    real(RP),         intent(inout) :: lst (IA,JA)
    real(RP),         intent(inout) :: albg(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP),         intent(inout) :: tc_urb(IA,JA)
    real(RP),         intent(inout) :: qc_urb(IA,JA)
    real(RP),         intent(inout) :: uc_urb(IA,JA)
    real(RP),         intent(inout) :: ust   (IA,JA)
    real(RP),         intent(inout) :: albu  (IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP),         intent(inout) :: lst_ocean(oIA_org,oJA_org)
    real(RP),         intent(inout) :: tw    (IA,JA)
    real(RP),         intent(inout) :: lz_org(lKA_org)
    real(RP),         intent(inout) :: topo_org(lIA_org,lJA_org)
    real(RP),         intent(inout) :: lmask_org(lIA_org,lJA_org)
    real(RP),         intent(inout) :: omask_org(oIA_org,oJA_org)
    real(RP),         intent(out)   :: sst (IA,JA)
    real(RP),         intent(out)   :: albw(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP),         intent(out)   :: z0w (IA,JA)
    character(len=*), intent(in)  :: basename_org
    character(len=*), intent(in)  :: basename_land
    character(len=*), intent(in)  :: basename_ocean
    integer,          intent(in)  :: mdlid_land
    integer,          intent(in)  :: mdlid_ocean
    integer,          intent(in)  :: ldims(3)
    integer,          intent(in)  :: odims(2)
    logical,          intent(in)  :: use_file_landwater                          ! use land water data from files
    real(RP),         intent(in)  :: init_landwater_ratio                        ! Ratio of land water to storage is constant,
!    real(RP),         intent(in)  :: init_landwater_ratio_each(LANDUSE_PFT_nmax) ! Ratio of land water to storage is constant,
                                                                                 ! if use_file_landwater is ".false." (each PFT)
    real(RP),         intent(in)  :: init_ocean_alb_lw
    real(RP),         intent(in)  :: init_ocean_alb_sw
    real(RP),         intent(in)  :: init_ocean_z0w
    integer,          intent(in)  :: intrp_iter_max
    logical,          intent(in)  :: soilwater_ds2vc_flag
    logical,          intent(in)  :: elevation_correction_land
    logical,          intent(in)  :: elevation_correction_ocean
    integer,          intent(in)  :: oistep
    integer,          intent(in)  :: listep
    logical,          intent(in)  :: multi_land
    logical,          intent(in)  :: URBAN_do

   ! land
    real(RP) :: tg_org   (lKA_org,lIA_org,lJA_org)
    real(RP) :: strg_org (lKA_org,lIA_org,lJA_org)
    real(RP) :: smds_org (lKA_org,lIA_org,lJA_org)
    real(RP) :: lst_org  (        lIA_org,lJA_org)
    real(RP) :: ust_org  (        lIA_org,lJA_org)
    real(RP) :: albg_org (        lIA_org,lJA_org,N_RAD_DIR,N_RAD_RGN)

    ! ocean
    real(RP) :: tw_org   (oIA_org,oJA_org)
    real(RP) :: sst_org  (oIA_org,oJA_org)
    real(RP) :: z0w_org  (oIA_org,oJA_org)
    real(RP) :: albw_org (oIA_org,oJA_org,N_RAD_DIR,N_RAD_RGN)
    real(RP) :: omask    (oIA_org,oJA_org)

    ! elevation correction
    real(RP) :: work(lIA_org,lJA_org)

    integer :: i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('___SurfaceInput',3)

    if ( do_read_land .and. ( first_surface .or. multi_land ) ) then

       select case( mdlid_land )
       case( iNetCDF ) ! TYPE: NetCDF

          call ParentLandInputNetCDF( &
               lKA_org, lKS_org, lKE_org, &
               lIA_org, lIS_org, lIE_org, &
               lJA_org, lJS_org, lJE_org, &
               tg_org, strg_org,           & ! (out)
               lst_org, ust_org, albg_org, & ! (out)
               topo_org, lmask_org,        & ! (inout)
               lz_org,                     & ! (inout)
               use_file_landwater,         & ! (in)
               ldims,                      & ! (in)
               listep                          ) ! (in)

       case( iGrADS ) ! TYPE: GrADS format

          call ParentLandInputGrADS( &
               lKA_org, lKS_org, lKE_org, &
               lIA_org, lIS_org, lIE_org, &
               lJA_org, lJS_org, lJE_org, &
               tg_org, strg_org, smds_org, & ! (out)
               lst_org,                    & ! (out)
               lz_org,                     & ! (out)
               topo_org, lmask_org,        & ! (out)
               use_waterratio,             & ! (out)
               ldims,                      & ! (in)
               basename_land,              & ! (in)
               use_file_landwater,         & ! (in)
               listep                      ) ! (in)
          ust_org = UNDEF
          albg_org = UNDEF

       end select

    end if

    call PROF_rapend  ('___SurfaceInput',3)

    call PROF_rapstart('___SurfaceBcast',3)

    if ( serial_land .and. ( first_surface .or. multi_land ) ) then
       call COMM_bcast( lKA_org, lIA_org, lJA_org, tg_org )
       if ( use_waterratio ) then
          call COMM_bcast( lKA_org, lIA_org, lJA_org, smds_org )
       else
          call COMM_bcast( lKA_org, lIA_org, lJA_org, strg_org )
       end if
       call COMM_bcast( lIA_org, lJA_org, lst_org )
       if ( URBAN_do ) call COMM_bcast( lIA_org, lJA_org, ust_org )
       call COMM_bcast( lIA_org, lJA_org, N_RAD_DIR, N_RAD_RGN, albg_org(:,:,:,:) )
       call COMM_bcast( lIA_org, lJA_org, topo_org )
       call COMM_bcast( lIA_org, lJA_org, lmask_org )
       call COMM_bcast( lKA_org, lz_org )
    end if

    call PROF_rapend  ('___SurfaceBcast',3)

    call PROF_rapstart('___SurfaceInput',3)

    if ( do_read_ocean ) then

       select case( mdlid_ocean )
       case( iNetCDF ) ! TYPE: NetCDF

          call ParentOceanInputNetCDF( &
               oIA_org, oIS_org, oIE_org, &
               oJA_org, oJS_org, oJE_org, &
               tw_org, sst_org,   & ! (out)
               albw_org, z0w_org, & ! (out)
               omask_org,         & ! (inout)
               odims,             & ! (in)
               oistep             ) ! (in)

       case( iGrADS ) ! TYPE: GrADS format

          call ParentOceanInputGrADS( &
               oIA_org, oIS_org, oIE_org, &
               oJA_org, oJS_org, oJE_org, &
               tw_org, sst_org,       & ! (out)
               omask_org,             & ! (out)
               basename_ocean, odims, & ! (in)
               oistep                 ) ! (in)
          albw_org = UNDEF
          z0w_org = UNDEF

       end select

    end if

    call PROF_rapend  ('___SurfaceInput',3)

    call PROF_rapstart('___SurfaceBcast',3)

    if ( serial_ocean ) then
       call COMM_bcast( oIA_org, oJA_org, tw_org )
       call COMM_bcast( oIA_org, oJA_org, sst_org )
       call COMM_bcast( oIA_org, oJA_org, N_RAD_DIR, N_RAD_RGN, albw_org(:,:,:,:) )

       call COMM_bcast( oIA_org, oJA_org, z0w_org )
       call COMM_bcast( oIA_org, oJA_org, omask_org )
    end if

    call PROF_rapend  ('___SurfaceBcast',3)

    call PROF_rapstart('___SurfaceInterp',3)

    if ( first_surface ) then

       if (    lIA_org .ne. oIA_org &
          .or. lJA_org .ne. oJA_org ) then
          ol_interp = .true.
       else
          ol_interp = .false.
   outer: do j = 1, lJA_org
          do i = 1, lIA_org
             if (    llon_org(i,j) .ne. olon_org(i,j) &
                .or. llat_org(i,j) .ne. olat_org(i,j) ) then
                ol_interp = .true.
                exit outer
             end if
          end do
          end do outer
       end if

       if ( ol_interp ) then
          ! interpolation factor between outer ocean grid and land grid
          call INTERP_factor2d( itp_nh_ol,        & ! [IN]
                                lIA_org, lJA_org, & ! [IN]
                                oIA_org, oJA_org, & ! [IN]
                                llon_org(:,:),    & ! [IN]
                                llat_org(:,:),    & ! [IN]
                                olon_org(:,:),    & ! [IN]
                                olat_org(:,:),    & ! [IN]
                                igrd_ol (:,:,:),  & ! [OUT]
                                jgrd_ol (:,:,:),  & ! [OUT]
                                hfact_ol(:,:,:)   ) ! [OUT]
       end if
    end if


    ! Ocean temp: interpolate over the land
    if ( i_INTRP_OCEAN_TEMP .ne. i_intrp_off ) then
       select case( i_INTRP_OCEAN_TEMP )
       case( i_intrp_mask )
          call make_mask( omask, tw_org, oIA_org, oJA_org, landdata=.false.)
          !$omp parallel do
          do j = 1, oJA_org
          do i = 1, oIA_org
             if ( omask_org(i,j) .ne. UNDEF ) omask(i,j) = omask_org(i,j)
          end do
          end do
       case( i_intrp_fill )
          call make_mask( omask, tw_org, oIA_org, oJA_org, landdata=.false.)
       end select
       call interp_OceanLand_data(tw_org, omask, oIA_org, oJA_org, .false., intrp_iter_max)
    end if

    ! SST: interpolate over the land
    if ( i_INTRP_OCEAN_SFC_TEMP .ne. i_intrp_off ) then
       select case( i_INTRP_OCEAN_SFC_TEMP )
       case( i_intrp_mask )
          call make_mask( omask, sst_org, oIA_org, oJA_org, landdata=.false.)
          !$omp parallel do
          do j = 1, oJA_org
          do i = 1, oIA_org
             if ( omask_org(i,j) .ne. UNDEF ) omask(i,j) = omask_org(i,j)
          end do
          end do
       case( i_intrp_fill )
          call make_mask( omask, sst_org, oIA_org, oJA_org, landdata=.false.)
       end select
       call interp_OceanLand_data(sst_org, omask, oIA_org, oJA_org, .false., intrp_iter_max)
    end if

    if ( first_surface .or. multi_land ) then

       call land_interporation( &
            lKA_org, lIA_org, lJA_org, &
            oIA_org, oJA_org,          &
            tg(:,:,:), strg(:,:,:),     & ! (out)
            lst(:,:), albg(:,:,:,:),    & ! (out)
            tg_org, strg_org, smds_org, & ! (inout)
            lst_org, albg_org,          & ! (inout)
            sst_org,                    & ! (in)
            lmask_org,                  & ! (in)
            lsmask_nest,                & ! (in)
            topo_org,                   & ! (in)
            lz_org, llon_org, llat_org, & ! (in)
            LCZ, CX, CY, LON, LAT,      & ! (in)
            maskval_tg, maskval_strg,   & ! (in)
            init_landwater_ratio,       & ! (in)
            !               init_landwater_ratio_each(:),    & ! (in)
            use_file_landwater,         & ! (in)
            use_waterratio,             & ! (in)
            soilwater_ds2vc_flag,       & ! (in)
            elevation_correction_land,  & ! (in)
            intrp_iter_max,             & ! (in)
            ol_interp                   ) ! (in)

       !$omp parallel do
       do j = 1, lJA_org
       do i = 1, lIA_org
          if ( topo_org(i,j) > UNDEF ) then ! ignore UNDEF value
             work(i,j) = lst_org(i,j) + topo_org(i,j) * LAPS
          else
             work(i,j) = lst_org(i,j)
          end if
       end do
       end do

       if ( ol_interp ) then
          ! land surface temperature at ocean grid
          call INTERP_interp2d( itp_nh_ol,          & ! [IN]
                                lIA_org, lJA_org, & ! [IN]
                                oIA_org, oJA_org, & ! [IN]
                                igrd_ol  (:,:,:),   & ! [IN]
                                jgrd_ol  (:,:,:),   & ! [IN]
                                hfact_ol (:,:,:),   & ! [IN]
                                work     (:,:),     & ! [IN]
                                lst_ocean(:,:)      ) ! [OUT]
       else
          !$omp parallel do
          do j = 1, oJA_org
          do i = 1, oIA_org
             lst_ocean(i,j) = work(i,j)
          end do
          end do
       end if

    end if

    call replace_misval_map( sst_org, lst_ocean, oIA_org, oJA_org, "SST" )
    call replace_misval_map( tw_org,  lst_ocean, oIA_org, oJA_org, "OCEAN_TEMP" )

    call ocean_interporation( oIA_org, oJA_org,                     & ! (in)
                              sst_org(:,:), tw_org(:,:),            & ! (in)
                              albw_org(:,:,:,:), z0w_org(:,:),      & ! (inout)
                              CX(:), CY(:),                         & ! (in)
                              elevation_correction_ocean,           & ! (in)
                              init_ocean_alb_lw, init_ocean_alb_sw, & ! (in)
                              init_ocean_z0w,                       & ! (in)
                              first_surface,                        & ! (in)
                              sst(:,:), tw(:,:),                    & ! (out)
                              albw(:,:,:,:), z0w(:,:)               ) ! (out)


    if ( first_surface .or. multi_land ) then
       ! replace values over the ocean ####
       !$omp parallel do
       do j = 1, JA
       do i = 1, IA
          if( abs(lsmask_nest(i,j)-0.0_RP) < EPS ) then ! ocean grid
             lst(i,j) = sst(i,j)
          endif
       enddo
       enddo
    end if


    if ( URBAN_do .and. first_surface ) then
       call urban_input( lst(:,:), albg(:,:,:,:),               & ! [IN]
                         tc_urb(:,:), qc_urb(:,:), uc_urb(:,:), & ! [OUT]
                         ust(:,:), albu(:,:,:,:)                ) ! [OUT]
    end if


    first_surface = .false.

    call PROF_rapend  ('___SurfaceInterp',3)


    return
  end subroutine ParentSurfaceInput

  !-----------------------------------------------------------------------------
  !> Boundary Data Write
  subroutine BoundarySurfaceOutput( &
       tg, &
       strg, &
       lst, &
       tw, &
       sst, &
       z0, &
       multi_ocean, &
       multi_land, &
       fid, &
       vid, &
       timeintv, &
       istep )
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

    real(RP),         intent(in) :: tg(LKA,IA,JA,1)
    real(RP),         intent(in) :: strg(LKA,IA,JA,1)
    real(RP),         intent(in) :: lst(IA,JA,1)
    real(RP),         intent(in) :: tw(1,IA,JA,1)
    real(RP),         intent(in) :: sst(IA,JA,1)
    real(RP),         intent(in) :: z0(IA,JA,1)
    logical,          intent(in) :: multi_ocean
    logical,          intent(in) :: multi_land
    integer,          intent(in) :: fid
    integer,          intent(in) :: vid(10)
    real(DP),         intent(in) :: timeintv
    integer,          intent(in) :: istep

    real(DP) :: timeofs
    !---------------------------------------------------------------------------

    call PROF_rapstart('___SurfaceOutput',3)

    timeofs = (istep - 1) * timeintv

    if ( multi_land ) then
       call FILE_CARTESC_write_var( fid, vid(1), tg(:,:,:,:),   'LAND_TEMP',     'LXYT', timeintv, timeofs=timeofs )
       call FILE_CARTESC_write_var( fid, vid(2), strg(:,:,:,:), 'LAND_WATER',    'LXYT', timeintv, timeofs=timeofs )
       call FILE_CARTESC_write_var( fid, vid(3), lst(:,:,:),    'LAND_SFC_TEMP', 'XYT',  timeintv, timeofs=timeofs )
    end if

    if ( multi_ocean ) then
       call FILE_CARTESC_write_var( fid, vid(6),  tw(:,:,:,:), 'OCEAN_TEMP',     'OXYT', timeintv, timeofs=timeofs )
       call FILE_CARTESC_write_var( fid, vid(7),  sst(:,:,:),  'OCEAN_SFC_TEMP', 'XYT',  timeintv, timeofs=timeofs )
       call FILE_CARTESC_write_var( fid, vid(10), z0(:,:,:), '  OCEAN_SFC_Z0',   'XYT',  timeintv, timeofs=timeofs )
    end if

    call PROF_rapend  ('___SurfaceOutput',3)

    return
  end subroutine BoundarySurfaceOutput


  !-------------------------------
  subroutine land_interporation( &
       kmax, imax, jmax, oimax,ojmax, &
       tg, strg, lst, albg,        &
       tg_org, strg_org, smds_org, &
       lst_org, albg_org,          &
       sst_org,                    &
       lmask_org, lsmask_nest,     &
       topo_org,                   &
       lz_org, llon_org, llat_org, &
       LCZ, CX, CY,                &
       LON, LAT,                   &
       maskval_tg, maskval_strg,   &
       init_landwater_ratio,       &
!       init_landwater_ratio_each,  &
       use_file_landwater,         &
       use_waterratio,             &
       soilwater_ds2vc_flag,       &
       elevation_correction,       &
       intrp_iter_max,             &
       ol_interp                   )
    use scale_prc, only: &
         PRC_abort
    use scale_const, only: &
         UNDEF => CONST_UNDEF, &
         EPS   => CONST_EPS,   &
         I_SW  => CONST_I_SW,  &
         I_LW  => CONST_I_LW,  &
         PI    => CONST_PI,    &
         LAPS  => CONST_LAPS
    use scale_interp, only: &
         INTERP_factor2d, &
         INTERP_factor3d, &
         INTERP_interp2d, &
         INTERP_interp3d
    use scale_mapprojection, only: &
         MAPPROJECTION_lonlat2xy
    use scale_comm_cartesC, only: &
         COMM_vars8, &
         COMM_wait
    use scale_filter, only: &
         FILTER_hyperdiff
    use scale_topography, only: &
         TOPOGRAPHY_Zsfc
    use scale_landuse, only: &
         LANDUSE_PFT_nmax,  &
         LANDUSE_index_PFT
    use mod_land_vars, only: &
         convert_WS2VWC
    implicit none
    integer,  intent(in)    :: kmax, imax, jmax
    integer,  intent(in)    :: oimax, ojmax
    real(RP), intent(out)   :: tg(LKMAX,IA,JA)
    real(RP), intent(out)   :: strg(LKMAX,IA,JA)
    real(RP), intent(out)   :: lst(IA,JA)
    real(RP), intent(out)   :: albg(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(inout) :: tg_org(kmax,imax,jmax)
    real(RP), intent(inout) :: strg_org(kmax,imax,jmax)
    real(RP), intent(inout) :: smds_org(kmax,imax,jmax)
    real(RP), intent(inout) :: lst_org(imax,jmax)
    real(RP), intent(inout) :: albg_org(imax,jmax,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(inout) :: sst_org(oimax,ojmax)
    real(RP), intent(in)    :: lmask_org(imax,jmax)
    real(RP), intent(in)    :: lsmask_nest(IA,JA)
    real(RP), intent(in)    :: topo_org(imax,jmax)
    real(RP), intent(in)    :: lz_org(kmax)
    real(RP), intent(in)    :: llon_org(imax,jmax)
    real(RP), intent(in)    :: llat_org(imax,jmax)
    real(RP), intent(in)    :: LCZ(LKMAX)
    real(RP), intent(in)    :: CX(IA)
    real(RP), intent(in)    :: CY(JA)
    real(RP), intent(in)    :: LON(IA,JA)
    real(RP), intent(in)    :: LAT(IA,JA)
    real(RP), intent(in)    :: maskval_tg
    real(RP), intent(in)    :: maskval_strg
    real(RP), intent(in)    :: init_landwater_ratio
!    real(RP), intent(in)    :: init_landwater_ratio_each(LANDUSE_PFT_nmax)
    logical,  intent(in)    :: use_file_landwater
    logical,  intent(in)    :: use_waterratio
    logical,  intent(in)    :: soilwater_ds2vc_flag
    logical,  intent(in)    :: elevation_correction
    integer,  intent(in)    :: intrp_iter_max
    logical,  intent(in)    :: ol_interp

    real(RP) :: lmask(imax,jmax)
    real(RP) :: smds(LKMAX,IA,JA)

    ! data for interporation
    real(RP) :: hfact_l(imax,jmax,itp_nh_ol)
    integer  :: igrd_l (imax,jmax,itp_nh_ol)
    integer  :: jgrd_l (imax,jmax,itp_nh_ol)
    real(RP) :: lX_org (imax,jmax)
    real(RP) :: lY_org (imax,jmax)
    logical  :: zonal, pole
    integer  :: igrd (         IA,JA,itp_nh_l)
    integer  :: jgrd (         IA,JA,itp_nh_l)
    real(RP) :: hfact(         IA,JA,itp_nh_l)
    integer  :: kgrdl (LKMAX,2,IA,JA,itp_nh_l)
    real(RP) :: vfactl(LKMAX,  IA,JA,itp_nh_l)


    real(RP) :: sst_land(imax,jmax)
    real(RP) :: work (imax,jmax)
    real(RP) :: work2(imax,jmax)

    real(RP) :: lz3d_org(kmax,imax,jmax)
    real(RP) :: lcz_3D(LKMAX,IA,JA)

    ! elevation correction
    real(RP) :: topo(IA,JA)
    real(RP) :: tdiff

    real(RP) :: one2d(IA,JA)
    real(RP) :: one3d(LKMAX,IA,JA)

    integer :: k, i, j, m, n


    ! Surface skin temp: interpolate over the ocean
    if ( i_INTRP_LAND_SFC_TEMP .ne. i_intrp_off ) then
       select case( i_INTRP_LAND_SFC_TEMP )
       case( i_intrp_mask )
          call make_mask( lmask, lst_org, imax, jmax, landdata=.true.)
          !$omp parallel do
          do j = 1, jmax
          do i = 1, imax
             if ( lmask_org(i,j) .ne. UNDEF ) lmask(i,j) = lmask_org(i,j)
          end do
          end do
       case( i_intrp_fill )
          call make_mask( lmask, lst_org, imax, jmax, landdata=.true.)
       case default
          LOG_ERROR("land_interporation",*) 'INTRP_LAND_SFC_TEMP is invalid.'
          call PRC_abort
       end select
       call interp_OceanLand_data(lst_org, lmask, imax, jmax, .true., intrp_iter_max)
    end if

    if ( ol_interp ) then
       ! interpolation facter between outer land grid and ocean grid
       call INTERP_factor2d( itp_nh_ol,       & ! [IN]
                             oimax, ojmax,    & ! [IN]
                             imax, jmax,      & ! [IN]
                             olon_org(:,:),   & ! [IN]
                             olat_org(:,:),   & ! [IN]
                             llon_org(:,:),   & ! [IN]
                             llat_org(:,:),   & ! [IN]
                             igrd_l  (:,:,:), & ! [OUT]
                             jgrd_l  (:,:,:), & ! [OUT]
                             hfact_l (:,:,:)  ) ! [OUT]

       ! sst on land grid
       call INTERP_interp2d( itp_nh_ol,        & ! [IN]
                             oimax, ojmax,     & ! [IN]
                             imax, jmax,       & ! [IN]
                             igrd_l   (:,:,:), & ! [IN]
                             jgrd_l   (:,:,:), & ! [IN]
                             hfact_l  (:,:,:), & ! [IN]
                             sst_org  (:,:),   & ! [IN]
                             sst_land (:,:)    ) ! [OUT]
    else
       !$omp parallel do
       do j = 1, jmax
       do i = 1, imax
          sst_land(i,j) = sst_org(i,j)
       end do
       end do
    end if

    !$omp parallel do
    do j = 1, jmax
    do i = 1, imax
       if ( topo_org(i,j) > UNDEF + EPS ) then ! ignore UNDEF value
          sst_land(i,j) = sst_land(i,j) - topo_org(i,j) * LAPS
       end if
    end do
    end do

    call replace_misval_map( lst_org, sst_land, imax, jmax, "SKINT" )

    ! replace missing value
    !$omp parallel do
    do j = 1, jmax
    do i = 1, imax
!       if ( skinw_org(i,j) == UNDEF ) skinw_org(i,j) = 0.0_RP
!       if ( snowq_org(i,j) == UNDEF ) snowq_org(i,j) = 0.0_RP
!       if ( snowt_org(i,j) == UNDEF ) snowt_org(i,j) = TEM00
       do m = 1, N_RAD_DIR
          if( albg_org(i,j,m,I_R_IR ) == UNDEF ) albg_org(i,j,m,I_R_IR ) = 0.03_RP ! emissivity of general ground surface : 0.95-0.98
          if( albg_org(i,j,m,I_R_NIR) == UNDEF ) albg_org(i,j,m,I_R_NIR) = 0.22_RP
          if( albg_org(i,j,m,I_R_VIS) == UNDEF ) albg_org(i,j,m,I_R_VIS) = 0.22_RP
       end do
    end do
    end do

    ! Land temp: interpolate over the ocean
    if ( i_INTRP_LAND_TEMP .ne. i_intrp_off ) then
       do k = 1, kmax
          !$omp parallel do
          do j = 1, jmax
          do i = 1, imax
             work(i,j) = tg_org(k,i,j)
          end do
          end do
          select case( i_INTRP_LAND_TEMP )
          case( i_intrp_mask )
             call make_mask( lmask, work, imax, jmax, landdata=.true.)
             !$omp parallel do
             do j = 1, jmax
             do i = 1, imax
                if ( lmask_org(i,j) .ne. UNDEF ) lmask(i,j) = lmask_org(i,j)
             end do
             end do
          case( i_intrp_fill )
             call make_mask( lmask, work, imax, jmax, landdata=.true.)
          end select
          call interp_OceanLand_data( work, lmask, imax, jmax, .true., intrp_iter_max )
          !replace land temp using skin temp
          call replace_misval_map( work, lst_org, imax,  jmax,  "STEMP")
          !$omp parallel do
          do j = 1, jmax
          do i = 1, imax
             tg_org(k,i,j) = work(i,j)
          end do
          end do
       end do
    end if


    ! fill grid data
    !$omp parallel do collapse(2)
    do j = 1, jmax
    do i = 1, imax
       lz3d_org(:,i,j) = lz_org(:)
    end do
    end do

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
       lcz_3D(:,i,j) = LCZ(:)
    enddo
    enddo

    select case( itp_type_l )
    case ( i_intrp_linear )

       if ( imax == 1 .or. jmax == 1 ) then
          LOG_ERROR("land_interporation",*) 'LINER interpolation requires nx, ny > 1'
          LOG_ERROR_CONT(*)                 'Use "DIST-WEIGHT" as INTRP_TYPE of PARAM_MKINIT_REAL_LAND'
          call PRC_abort
       end if

       !$omp parallel do
       do j = 1, jmax
       do i = 1, imax
          work(i,j) = sign( min( abs(llat_org(i,j)), PI * 0.499999_RP ), llat_org(i,j) )
       end do
       end do

       call MAPPROJECTION_lonlat2xy( imax, 1, imax, jmax, 1, jmax, &
                                     llon_org(:,:), work(:,:), & ! [IN]
                                     lX_org(:,:), lY_org(:,:)  ) ! [OUT]

       zonal = ( maxval(llon_org) - minval(llon_org) ) > 2.0_RP * PI * 0.9_RP
       pole = ( maxval(llat_org) > PI * 0.5_RP * 0.9_RP ) .or. ( minval(llat_org) < - PI * 0.5_RP * 0.9_RP )
       call INTERP_factor3d( kmax, 1, kmax,        & ! [IN]
                             imax, jmax,           & ! [IN]
                             LKMAX, LKS, LKE,      & ! [IN]
                             IA, JA,               & ! [IN]
                             lX_org(:,:),          & ! [IN]
                             lY_org(:,:),          & ! [IN]
                             lz3d_org(:,:,:),      & ! [IN]
                             CX(:), CY(:),         & ! [IN]
                             lcz_3D  (:,:,:),      & ! [IN]
                             igrd    (    :,:,:),  & ! [OUT]
                             jgrd    (    :,:,:),  & ! [OUT]
                             hfact   (    :,:,:),  & ! [OUT]
                             kgrdl   (:,:,:,:,:),  & ! [OUT]
                             vfactl  (:,  :,:,:),  & ! [OUT]
                             flag_extrap = .true., & ! [IN, optional]
                             zonal = zonal,        & ! [IN, optional]
                             pole  = pole          ) ! [IN, optional]

    case ( I_intrp_dstwgt )

       call INTERP_factor3d( itp_nh_l,            & ! [IN]
                             kmax, 1, kmax,       & ! [IN]
                             imax, jmax,          & ! [IN]
                             LKMAX, LKS, LKE,     & ! [IN]
                             IA, JA,              & ! [IN]
                             llon_org(:,:),       & ! [IN]
                             llat_org(:,:),       & ! [IN]
                             lz3d_org(:,:,:),     & ! [IN]
                             LON(:,:), LAT(:,:),  & ! [IN]
                             lcz_3D  (:,:,:),     & ! [IN]
                             igrd    (    :,:,:), & ! [OUT]
                             jgrd    (    :,:,:), & ! [OUT]
                             hfact   (    :,:,:), & ! [OUT]
                             kgrdl   (:,:,:,:,:), & ! [OUT]
                             vfactl  (:,  :,:,:), & ! [OUT]
                             flag_extrap = .true. ) ! [IN, optional]

    end select

    call INTERP_interp2d( itp_nh_l,        & ! [IN]
                          imax, jmax,      & ! [IN]
                          IA, JA,          & ! [IN]
                          igrd    (:,:,:), & ! [IN]
                          jgrd    (:,:,:), & ! [IN]
                          hfact   (:,:,:), & ! [IN]
                          lst_org (:,:),   & ! [IN]
                          lst     (:,:)    ) ! [OUT]
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
                              lst(:,:), FILTER_ORDER, FILTER_NITER )
       call COMM_vars8( lst(:,:), 1 )
       call COMM_wait ( lst(:,:), 1, .false. )
    end if


    if ( FILTER_NITER > 0 ) then
       !$omp parallel do
       do j = 1, JA
       do i = 1, IA
          one2d(i,j) = 1.0_RP
       end do
       end do
    end if

    do n = 1, N_RAD_RGN
    do m = 1, N_RAD_DIR

       call INTERP_interp2d( itp_nh_l,          & ! [IN]
                             imax, jmax,        & ! [IN]
                             IA, JA,            & ! [IN]
                             igrd    (:,:,:),   & ! [IN]
                             jgrd    (:,:,:),   & ! [IN]
                             hfact   (:,:,:),   & ! [IN]
                             albg_org(:,:,m,n), & ! [IN]
                             albg    (:,:,m,n)  ) ! [OUT]
       if ( FILTER_NITER > 0 ) then
          call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
                                 albg(:,:,m,n), FILTER_ORDER, FILTER_NITER, &
                                 limiter_sign = one2d(:,:) )
          call COMM_vars8( albg(:,:,m,n), 1 )
          call COMM_wait ( albg(:,:,m,n), 1, .false. )
       end if
    end do
    end do

    call INTERP_interp3d( itp_nh_l,            &
                          kmax, 1, kmax,       &
                          imax, jmax,          &
                          LKMAX, LKS, LKE,     &
                          IA, JA,              &
                          igrd    (    :,:,:), & ! [IN]
                          jgrd    (    :,:,:), & ! [IN]
                          hfact   (    :,:,:), & ! [IN]
                          kgrdl   (:,:,:,:,:), & ! [IN]
                          vfactl  (:,  :,:,:), & ! [IN]
                          lz3d_org(:,:,:),     & ! [IN]
                          lcz_3D  (:,:,:),     & ! [IN]
                          tg_org  (:,:,:),     & ! [IN]
                          tg      (:,:,:)      ) ! [OUT]

    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       tg(LKMAX,i,j) = tg(LKMAX-1,i,j)
    enddo
    enddo

    ! replace values over the ocean
    do k = 1, LKMAX
       call replace_misval_const( tg(k,:,:), maskval_tg, lsmask_nest )
    enddo
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( LKMAX, 1, LKMAX, IA, ISB, IEB, JA, JSB, JEB, &
                              tg(:,:,:), FILTER_ORDER, FILTER_NITER )
       call COMM_vars8( tg(:,:,:), 1 )
       call COMM_wait ( tg(:,:,:), 1, .false. )
    end if


    ! elevation correction
    if ( elevation_correction ) then
       call INTERP_interp2d( itp_nh_l,        & ! [IN]
                             imax, jmax,      & ! [IN]
                             IA, JA,          & ! [IN]
                             igrd    (:,:,:), & ! [IN]
                             jgrd    (:,:,:), & ! [IN]
                             hfact   (:,:,:), & ! [IN]
                             topo_org(:,:),   & ! [IN]
                             topo    (:,:)    ) ! [OUT]
       if ( FILTER_NITER > 0 ) then
          call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
                                 topo(:,:), FILTER_ORDER, FILTER_NITER )
          call COMM_vars8( topo(:,:), 1 )
          call COMM_wait ( topo(:,:), 1, .false. )
       end if

       !$omp parallel do &
       !$omp private(tdiff)
       do j = 1, JA
       do i = 1, IA
          if ( topo(i,j) > UNDEF + EPS ) then ! ignore UNDEF value
             tdiff = ( TOPOGRAPHY_Zsfc(i,j) - topo(i,j) ) * LAPS
             lst(i,j) = lst(i,j) - tdiff
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
             do k = 1, kmax
                !$omp parallel do
                do j = 1, jmax
                do i = 1, imax
                   work(i,j) = smds_org(k,i,j)
                end do
                end do
                select case( i_INTRP_LAND_WATER )
                case( i_intrp_mask )
                   call make_mask( lmask, work, imax, jmax, landdata=.true.)
                   !$omp parallel do
                   do j = 1, jmax
                   do i = 1, imax
                      if ( lmask_org(i,j) .ne. UNDEF ) lmask(i,j) = lmask_org(i,j)
                   end do
                   end do
                case( i_intrp_fill )
                   call make_mask( lmask, work, imax, jmax, landdata=.true.)
                end select
                call interp_OceanLand_data(work, lmask, imax, jmax, .true., intrp_iter_max)
                !$omp parallel do
                do j = 1, jmax
                do i = 1, imax
!                   work2(i,j) = init_landwater_ratio_each( LANDUSE_index_PFT(i,j,1) )
                   work2(i,j) = init_landwater_ratio
                end do
                end do
                !replace missing value to init_landwater_ratio_each
                call replace_misval_map( work, work2, imax, jmax, "SMOISDS")
                !$omp parallel do
                do j = 1, jmax
                do i = 1, imax
                   smds_org(k,i,j) = work(i,j)
                end do
                end do
             enddo
          end if

          call INTERP_interp3d( itp_nh_l,            &
                                kmax, 1, kmax,       &
                                imax, jmax,          &
                                LKMAX, LKS, LKE,     &
                                IA, JA,              &
                                igrd    (    :,:,:), & ! [IN]
                                jgrd    (    :,:,:), & ! [IN]
                                hfact   (    :,:,:), & ! [IN]
                                kgrdl   (:,:,:,:,:), & ! [IN]
                                vfactl  (:,  :,:,:), & ! [IN]
                                lz3d_org(:,:,:),     & ! [IN]
                                lcz_3D  (:,:,:),     & ! [IN]
                                smds_org(:,:,:),     & ! [IN]
                                smds    (:,:,:)      ) ! [OUT]

          do k = 1, LKMAX-1
             strg(k,:,:) = convert_WS2VWC( smds(k,:,:), critical=soilwater_DS2VC_flag )
          end do

       else

          if ( i_INTRP_LAND_WATER .ne. i_intrp_off ) then
             do k = 1, kmax
                !$omp parallel do
                do j = 1, jmax
                do i = 1, imax
                   work(i,j) = strg_org(k,i,j)
                end do
                end do
                select case( i_INTRP_LAND_WATER )
                case( i_intrp_mask )
                   call make_mask( lmask, work, imax, jmax, landdata=.true.)
                   !$omp parallel do
                   do j = 1, jmax
                   do i = 1, imax
                      if ( lmask_org(i,j) .ne. UNDEF ) lmask(i,j) = lmask_org(i,j)
                   end do
                   end do
                case( i_intrp_fill )
                   call make_mask( lmask, work, imax, jmax, landdata=.true.)
                end select
                call interp_OceanLand_data(work, lmask, imax, jmax, .true., intrp_iter_max)
                !$omp parallel do
                do j = 1, jmax
                do i = 1, imax
                   lmask(i,j) = maskval_strg
                end do
                end do
                !replace missing value to init_landwater_ratio
                call replace_misval_map( work, lmask, imax, jmax, "SMOIS")
                !$omp parallel do
                do j = 1, jmax
                do i = 1, imax
                   strg_org(k,i,j) = work(i,j)
                end do
                end do
             enddo
          end if

          call INTERP_interp3d( itp_nh_l,            &
                                kmax, 1, kmax,       &
                                imax, jmax,          &
                                LKMAX, LKS, LKE,     &
                                IA, JA,              &
                                igrd    (    :,:,:), & ! [IN]
                                jgrd    (    :,:,:), & ! [IN]
                                hfact   (    :,:,:), & ! [IN]
                                kgrdl   (:,:,:,:,:), & ! [IN]
                                vfactl  (:,  :,:,:), & ! [IN]
                                lz3d_org(:,:,:),     & ! [IN]
                                lcz_3D  (:,:,:),     & ! [IN]
                                strg_org(:,:,:),     & ! [IN]
                                strg    (:,:,:)      ) ! [OUT]
       end if

       ! replace values over the ocean
       do k = 1, LKMAX-1
          call replace_misval_const( strg(k,:,:), maskval_strg, lsmask_nest )
       enddo

       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = 1, LKMAX
          strg(k,i,j) = max( min( strg(k,i,j), 1.0_RP ), 0.0_RP )
       end do
       end do
       end do

       if ( FILTER_NITER > 0 ) then
          !$omp parallel do collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = 1, LKMAX-1
             one3d(k,i,j) = 1.0_RP
          end do
          end do
          end do
          call FILTER_hyperdiff( LKMAX, 1, LKMAX-1, IA, ISB, IEB, JA, JSB, JEB, &
                                 strg(:,:,:), FILTER_ORDER, FILTER_NITER, &
                                 limiter_sign = one3d(:,:,:) )
          call COMM_vars8( strg(:,:,:), 1 )
          call COMM_wait ( strg(:,:,:), 1, .false. )
       end if

       !$omp parallel do
       do j = 1, JA
       do i = 1, IA
          strg(LKMAX,i,j) = strg(LKMAX-1,i,j)
       enddo
       enddo

    else  ! not read from boundary file

       do k = 1, LKMAX
          !$omp parallel do
          do j = 1, JA
          do i = 1, IA
!             work(i,j) = init_landwater_ratio_each( LANDUSE_index_PFT(i,j,1) )
             work(i,j) = init_landwater_ratio
          end do
          end do
          ! conversion from water saturation [fraction] to volumetric water content [m3/m3]
          strg(k,:,:) = convert_WS2VWC( work(:,:), critical=soilwater_DS2VC_flag )
       end do

    endif ! use_file_waterratio


    return
  end subroutine land_interporation

  subroutine ocean_interporation( &
       imax, jmax, &
       sst_org, tw_org, albw_org, z0w_org,   &
       CX, CY,                               &
       elevation_correction_ocean,           &
       init_ocean_alb_lw, init_ocean_alb_sw, &
       init_ocean_z0w,                       &
       first_surface,                        &
       sst, tw, albw, z0w                    )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       PI    => CONST_PI,    &
       LAPS  => CONST_LAPS
    use scale_topography, only: &
       TOPOGRAPHY_Zsfc
    use scale_interp, only: &
       INTERP_factor2d, &
       INTERP_interp2d
    use scale_filter, only: &
       FILTER_hyperdiff
    use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    integer,  intent(in)    :: imax, jmax
    real(RP), intent(in)    :: sst_org (imax,jmax)
    real(RP), intent(in)    :: tw_org  (imax,jmax)
    real(RP), intent(inout) :: albw_org(imax,jmax,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(inout) :: z0w_org (imax,jmax)
    real(RP), intent(in)    :: CX(IA)
    real(RP), intent(in)    :: CY(JA)
    logical,  intent(in)    :: elevation_correction_ocean
    real(RP), intent(in)    :: init_ocean_alb_lw
    real(RP), intent(in)    :: init_ocean_alb_sw
    real(RP), intent(in)    :: init_ocean_z0w
    logical,  intent(in)    :: first_surface

    real(RP), intent(out) :: sst (IA,JA)
    real(RP), intent(out) :: tw  (IA,JA)
    real(RP), intent(out) :: albw(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(out) :: z0w (IA,JA)

    ! for interpolation
    real(RP) :: oX_org(imax,jmax)
    real(RP) :: oY_org(imax,jmax)
    logical  :: zonal, pole

    real(RP) :: one(IA,JA)
    real(RP) :: tdiff

    integer :: i, j, m, n

    !$omp parallel do
    do j = 1, jmax
    do i = 1, imax
       do m = 1, N_RAD_DIR
          if ( albw_org(i,j,m,I_R_IR ) == UNDEF ) albw_org(i,j,m,I_R_IR ) = init_ocean_alb_lw
          if ( albw_org(i,j,m,I_R_NIR) == UNDEF ) albw_org(i,j,m,I_R_NIR) = init_ocean_alb_sw
          if ( albw_org(i,j,m,I_R_VIS) == UNDEF ) albw_org(i,j,m,I_R_VIS) = init_ocean_alb_sw
          if ( albw_org(i,j,m,I_R_VIS) == UNDEF ) albw_org(i,j,m,I_R_VIS) = init_ocean_alb_sw
       end do
       if ( z0w_org(i,j) == UNDEF ) z0w_org(i,j) = init_ocean_z0w
    end do
    end do

    if ( first_surface ) then

       ! interporation for ocean variables

       select case( itp_type_a )
       case ( i_intrp_linear )

          if ( imax == 1 .or. jmax == 1 ) then
             LOG_ERROR("ocean_interporation",*) 'LINER interpolation requires nx, ny > 1'
             LOG_ERROR_CONT(*)                  'Use "DIST-WEIGHT" as INTRP_TYPE of PARAM_MKINIT_REAL_OCEAN'
             call PRC_abort
          end if

          !$omp parallel do
          do j = 1, jmax
          do i = 1, imax
             olat_org(i,j) = sign( min( abs(olat_org(i,j)), PI * 0.499999_RP ), olat_org(i,j) )
          end do
          end do

          call MAPPROJECTION_lonlat2xy( imax, 1, imax, jmax, 1, jmax, &
                                        olon_org(:,:), olat_org(:,:), & ! [IN]
                                        oX_org(:,:), oY_org(:,:)      ) ! [OUT]

          zonal = ( maxval(olon_org) - minval(olon_org) ) > 2.0_RP * PI * 0.9_RP
          pole = ( maxval(olat_org) > PI * 0.5_RP * 0.9_RP ) .or. ( minval(olat_org) < - PI * 0.5_RP * 0.9_RP )
          call INTERP_factor2d( imax, jmax,    & ! [IN]
                                IA, JA,        & ! [IN]
                                oX_org(:,:),   & ! [IN]
                                oY_org(:,:),   & ! [IN]
                                CX(:), CY(:),  & ! [IN]
                                oigrd (:,:,:), & ! [OUT]
                                ojgrd (:,:,:), & ! [OUT]
                                ohfact(:,:,:), & ! [OUT]
                                zonal = zonal, & ! [IN]
                                pole  = pole   ) ! [IN]

       case ( I_intrp_dstwgt )

          call INTERP_factor2d( itp_nh_o,           & ! [IN]
                                imax, jmax,         & ! [IN]
                                IA, JA,             & ! [IN]
                                olon_org(:,:),      & ! [IN]
                                olat_org(:,:),      & ! [IN]
                                LON(:,:), LAT(:,:), & ! [IN]
                                oigrd (:,:,:),      & ! [OUT]
                                ojgrd (:,:,:),      & ! [OUT]
                                ohfact(:,:,:)       ) ! [OUT]

       end select

    end if


    call INTERP_interp2d( itp_nh_o,      & ! [IN]
                          imax, jmax,    & ! [IN]
                          IA, JA,        & ! [IN]
                          oigrd (:,:,:), & ! [IN]
                          ojgrd (:,:,:), & ! [IN]
                          ohfact(:,:,:), & ! [IN]
                          tw_org(:,:),   & ! [IN]
                          tw    (:,:)    ) ! [OUT]
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
                              tw(:,:), FILTER_ORDER, FILTER_NITER )
       call COMM_vars8( tw(:,:), 1 )
       call COMM_wait ( tw(:,:), 1, .false. )
    end if

    call INTERP_interp2d( itp_nh_o,       & ! [IN]
                          imax, jmax,     & ! [IN]
                          IA, JA,         & ! [IN]
                          oigrd  (:,:,:), & ! [IN]
                          ojgrd  (:,:,:), & ! [IN]
                          ohfact (:,:,:), & ! [IN]
                          sst_org(:,:),   & ! [IN]
                          sst    (:,:)    ) ! [OUT]
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
                              sst(:,:), FILTER_ORDER, FILTER_NITER )
       call COMM_vars8( sst(:,:), 1 )
       call COMM_wait ( sst(:,:), 1, .false. )
    end if

    ! elevation correction
    if ( elevation_correction_ocean ) then

       !$omp parallel do &
       !$omp private(tdiff)
       do j = 1, JA
       do i = 1, IA
          tdiff = TOPOGRAPHY_Zsfc(i,j) * LAPS
          sst(i,j) = sst(i,j) - tdiff
          tw (i,j) = tw (i,j) - tdiff
       end do
       end do

    end if


    if ( FILTER_NITER > 0 ) then
       !$omp parallel do
       do j = 1, JA
       do i = 1, IA
          one(i,j) = 1.0_RP
       end do
       end do
    end if

    do n = 1, N_RAD_RGN
    do m = 1, N_RAD_DIR

       call INTERP_interp2d( itp_nh_o,          & ! [IN]
                             imax, jmax,        & ! [IN]
                             IA, JA,            & ! [IN]
                             oigrd   (:,:,:),   & ! [IN]
                             ojgrd   (:,:,:),   & ! [IN]
                             ohfact  (:,:,:),   & ! [IN]
                             albw_org(:,:,m,n), & ! [IN]
                             albw    (:,:,m,n)  ) ! [OUT]
       if ( FILTER_NITER > 0 ) then
          call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
                                 albw(:,:,m,n), FILTER_ORDER, FILTER_NITER, &
                                 limiter_sign = one(:,:) )
          call COMM_vars8( albw(:,:,m,n), 1 )
          call COMM_wait ( albw(:,:,m,n), 1, .false. )
       end if

    end do
    end do

    call INTERP_interp2d( itp_nh_o,        & ! [IN]
                          imax, jmax,      & ! [IN]
                          IA, JA,          & ! [IN]
                          oigrd   (:,:,:), & ! [IN]
                          ojgrd   (:,:,:), & ! [IN]
                          ohfact  (:,:,:), & ! [IN]
                          z0w_org (:,:),   & ! [IN]
                          z0w     (:,:)    ) ! [OUT]
    if ( FILTER_NITER > 0 ) then
       call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
                              z0w(:,:), FILTER_ORDER, FILTER_NITER, &
                              limiter_sign = one(:,:) )
       call COMM_vars8( z0w(:,:), 1 )
       call COMM_wait ( z0w(:,:), 1, .false. )
    end if


    return
  end subroutine ocean_interporation

  !-------------------------------
  subroutine urban_input( &
       lst, albg,              &
       tc_urb, qc_urb, uc_urb, &
       ust, albu               )
    use mod_atmos_vars, only: &
         DENS, &
         MOMX, &
         MOMY, &
         RHOT, &
         QTRC
    use scale_atmos_hydrometeor, only: &
         I_QV
    use scale_atmos_thermodyn, only: &
         THERMODYN_specific_heat  => ATMOS_THERMODYN_specific_heat, &
         THERMODYN_rhot2temp_pres => ATMOS_THERMODYN_rhot2temp_pres
    use scale_comm_cartesC, only: &
         COMM_vars8, &
         COMM_wait
    implicit none
    real(RP), intent(in)  :: lst   (IA,JA)
    real(RP), intent(in)  :: albg  (IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(out) :: tc_urb(IA,JA)
    real(RP), intent(out) :: qc_urb(IA,JA)
    real(RP), intent(out) :: uc_urb(IA,JA)
    real(RP), intent(out) :: ust   (IA,JA)
    real(RP), intent(out) :: albu  (IA,JA,N_RAD_DIR,N_RAD_RGN)

    real(RP) :: temp, pres
    real(RP) :: Qdry
    real(RP) :: Rtot
    real(RP) :: CVtot
    real(RP) :: CPtot

    integer :: i, j

    ! urban data

    !$omp parallel do collapse(2) &
    !$omp private(Qdry,Rtot,CVtot,CPtot,temp,pres)
    do j = 1, JA
    do i = 1, IA
       call THERMODYN_specific_heat( QA, &
                                     qtrc(KS,i,j,:), &
                                     TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! [IN]
                                     Qdry, Rtot, CVtot, CPtot                                 ) ! [OUT]
       call THERMODYN_rhot2temp_pres( dens(KS,i,j), rhot(KS,i,j), Rtot, CVtot, CPtot, & ! [IN]
                                      temp, pres                                      ) ! [OUT]

       tc_urb(i,j) = temp
       if ( I_QV > 0 ) then
          qc_urb(i,j) = qtrc(KS,i,j,I_QV)
       else
          qc_urb(i,j) = 0.0_RP
       end if
    enddo
    enddo

    !$omp parallel do
    do j = 1, JA-1
    do i = 1, IA-1
       uc_urb(i,j) = max(sqrt( ( momx(KS,i,j) / (dens(KS,i+1,  j)+dens(KS,i,j)) * 2.0_RP )**2.0_RP &
                             + ( momy(KS,i,j) / (dens(KS,  i,j+1)+dens(KS,i,j)) * 2.0_RP )**2.0_RP ), &
                             0.01_RP)
    enddo
    enddo
    !$omp parallel do
    do j = 1, JA-1
       uc_urb(IA,j) = max(sqrt( ( momx(KS,IA,j) /  dens(KS,IA,j  ) )**2.0_RP &
                              + ( momy(KS,IA,j) / (dens(KS,IA,j+1)+dens(KS,IA,j)) * 2.0_RP )**2.0_RP ), &
                              0.01_RP)
    enddo
    !$omp parallel do
    do i = 1, IA-1
       uc_urb(i,JA) = max(sqrt( ( momx(KS,i,JA) / (dens(KS,i+1,JA)+dens(KS,i,JA)) * 2.0_RP )**2.0_RP &
                              + ( momy(KS,i,JA) /  dens(KS,i  ,JA) )**2.0_RP ), 0.01_RP)
    enddo
    uc_urb(IA,JA) = max(sqrt( ( momx(KS,IA,JA) / dens(KS,IA,JA) )**2.0_RP &
                            + ( momy(KS,IA,JA) / dens(KS,IA,JA) )**2.0_RP ), 0.01_RP)

    call COMM_vars8( uc_urb, 1 )
    call COMM_wait ( uc_urb, 1, .false. )


!!$    ! Urban surface temp: interpolate over the ocean
!!$    if ( i_INTRP_URB_SFC_TEMP .ne. i_intrp_off ) then
!!$       select case( i_INTRP_URB_SFC_TEMP )
!!$       case( i_intrp_mask )
!!$          call make_mask( lmask, ust_org, imax, jmax, landdata=.true.)
!!$          !$omp parallel do
!!$          do j = 1, jmax
!!$          do i = 1, imax
!!$             if ( lmask_org(i,j) .ne. UNDEF ) lmask(i,j) = lmask_org(i,j)
!!$          end do
!!$          end do
!!$       case( i_intrp_fill )
!!$          call make_mask( lmask, ust_org, imax, jmax, landdata=.true.)
!!$       case default
!!$          LOG_ERROR("urban_input",*) 'INTRP_URB_SFC_TEMP is invalid.'
!!$          call PRC_abort
!!$       end select
!!$       call interp_OceanLand_data(ust_org, lmask, imax, jmax, .true., intrp_iter_max)
!!$    end if
!!$
!!$    !$omp parallel do
!!$    do j = 1, jmax
!!$    do i = 1, imax
!!$       if ( ust_org(i,j) == UNDEF ) ust_org(i,j) = lst_org(i,j)
!!$    end do
!!$    end do
!!$
!!$    call INTERP_interp2d( itp_nh_l,           & ! [IN]
!!$                          imax, jmax, & ! [IN]
!!$                          IA, JA,             & ! [IN]
!!$                          igrd    (:,:,:),    & ! [IN]
!!$                          jgrd    (:,:,:),    & ! [IN]
!!$                          hfact   (:,:,:),    & ! [IN]
!!$                          ust_org (:,:),      & ! [IN]
!!$                          ust     (:,:)       ) ! [OUT]
!!$    if ( FILTER_NITER > 0 ) then
!!$       call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
!!$                              ust(:,:), FILTER_ORDER, FILTER_NITER )
!!$       call COMM_vars8( ust(:,:), 1 )
!!$       call COMM_wait ( ust(:,:), 1, .false. )
!!$    end if
!!$
!!$    !$omp parallel do
!!$    do j = 1, JA
!!$    do i = 1, IA
!!$       if( abs(lsmask_nest(i,j)-0.0_RP) < EPS ) then ! ocean grid
!!$          ust(i,j) = sst(i,j,nn)
!!$       endif
!!$    enddo
!!$    enddo
!!$
!!$    !$omp parallel do &
!!$    !$omp private(tdiff)
!!$    do j = 1, JA
!!$    do i = 1, IA
!!$       if ( topo(i,j) > 0.0_RP ) then ! ignore UNDEF value
!!$          tdiff = ( TOPOGRAPHY_Zsfc(i,j) - topo(i,j) ) * LAPS
!!$          ust(i,j) = ust(i,j) - tdiff
!!$       end if
!!$    end do
!!$    end do


    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       ust(i,j) = lst(i,j)
    end do
    end do


    ! copy albedo of land to urban
    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       albu(i,j,:,:) = albg(i,j,:,:)
    enddo
    enddo

    return
  end subroutine urban_input

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
       !$omp parallel do
       do j = 1, ny
       do i = 1, nx
          gmask(i,j) = 1.0_RP  ! gmask=1 will be skip in "interp_OceanLand_data"
       end do
       end do
       dd         = 0.0_RP
    else
       !$omp parallel do
       do j = 1, ny
       do i = 1, nx
          gmask(i,j) = 0.0_RP  ! gmask=0 will be skip in "interp_OceanLand_data"
       end do
       end do
       dd         = 1.0_RP
    endif

    !$omp parallel do
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
       UNDEF => CONST_UNDEF, &
       EPS   => CONST_EPS
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

    LOG_NEWLINE
    LOG_INFO("interp_OceanLand_data",*) 'Interpolation'

    if ( landdata ) then
       LOG_INFO("interp_OceanLand_data",*) 'target mask : LAND'
       mask_target = 1 ! interpolation for land data
    else
       LOG_INFO("interp_OceanLand_data",*) 'target mask : OCEAN'
       mask_target = 0 ! interpolation for ocean data
    endif

    ! search target cell for interpolation
    num_land  = 0
    num_ocean = 0
    !$omp parallel do &
    !$omp reduction(+:num_land,num_ocean)
    do j = 1, ny
    do i = 1, nx
       mask(i,j) = int( 0.5_RP - sign(0.5_RP,abs(lsmask(i,j)-1.0_RP)-EPS) ) ! 1 for land, 0 for ocean
       num_land  = num_land  + (   mask(i,j) )
       num_ocean = num_ocean + ( 1-mask(i,j) )
    enddo
    enddo

    LOG_PROGRESS('(1x,A,I3.3,A,2I8)') 'ite=', 0, ', (land,ocean) = ', num_land, num_ocean

    ! start interpolation
    do ite = 1, iter_max
       ! save previous state
       !$omp parallel do
       do j = 1, ny
       do i = 1, nx
          mask_prev(i,j) = mask(i,j)
          data_prev(i,j) = data(i,j)
       end do
       end do
       num_replaced = 0

       !$omp parallel do collapse(2) &
       !$omp private(istr,iend,jstr,jend,tmp,cnt,sw) &
       !$omp reduction(+:num_replaced)
       do j = 1, ny
       do i = 1, nx

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
!       LOG_PROGRESS('(1x,A,I3.3,A,3I8,A,I8)') 'ite=', ite, &
!                                              ', (land,ocean,replaced) = ', num_land, num_ocean, num_replaced, ' / ', nx*ny

       if( num_replaced == 0 ) exit

    enddo ! itelation

    LOG_PROGRESS('(1x,A,I3.3,A,2I8)') 'ite=', ite, ', (land,ocean) = ', num_land, num_ocean

    !$omp parallel do
    do j = 1, ny
    do i = 1, nx
       if ( abs(mask(i,j)-mask_target) > EPS ) data(i,j) = UNDEF
    end do
    end do


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

    !$omp parallel do
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
    logical :: error

    error = .false.
    !$omp parallel do
    do j = 1, ny
       if ( error ) cycle
       do i = 1, nx
          if( abs(data(i,j) - UNDEF) < sqrt(EPS) )then
             if( abs(maskval(i,j) - UNDEF) < sqrt(EPS) )then
                LOG_ERROR("replace_misval_map",*) "data for mask of "//trim(elem)//"(",i,",",j,") includes missing value."
                error = .true.
                exit
             else
                data(i,j) = maskval(i,j)
             endif
          endif
       enddo
    enddo

    if ( error ) then
       LOG_ERROR_CONT(*) "Please check input data of SKINTEMP or SST. "
       call PRC_abort
    end if

    return
  end subroutine replace_misval_map

  subroutine get_IJrange( &
       IS_org, IE_org, JS_org, JE_org,      &
       IA_org, JA_org,                     &
       LON_min, LON_max, LAT_min, LAT_max, &
       LON_all, LAT_all                    )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LON, &
       ATMOS_GRID_CARTESC_REAL_LAT
    integer, intent(out) :: IS_org
    integer, intent(out) :: IE_org
    integer, intent(out) :: JS_org
    integer, intent(out) :: JE_org

    integer,  intent(in) :: IA_org
    integer,  intent(in) :: JA_org
    real(RP), intent(in) :: LON_min, LON_max
    real(RP), intent(in) :: LAT_min, LAT_max
    real(RP), intent(in) :: LON_all(IA_org,JA_org)
    real(RP), intent(in) :: LAT_all(IA_org,JA_org)

    real(RP) :: min, max

    logical :: LON_mask(IA_org)
    logical :: LAT_mask(JA_org)

    integer :: i, j

    if ( LON_min < minval( LON_all ) .or. LON_max > maxval( LON_all ) ) then
       ! probably global (cyclic) data
       IS_org = 1
       IE_org = IA_org
    else
       min = maxval( minval( LON_all(:,:), dim=2 ), mask=all( LON_all(:,:) < LON_min, dim=2 ) )
       max = minval( maxval( LON_all(:,:), dim=2 ), mask=all( LON_all(:,:) > LON_max, dim=2 ) )
       LON_mask(:) = any( LON_all(:,:) - min > -EPS, dim=2 ) .AND. any( LON_all(:,:) - max < EPS, dim=2 )
       do i = 1, IA_org
          if( LON_mask(i) ) then
             IS_org = i
             exit
          endif
       end do
       do i = IA_org, 1, -1
          if( LON_mask(i) ) then
             IE_org = i
             exit
          endif
       end do
    end if

    if ( LAT_min < minval( LAT_all ) .or. LAT_max > maxval( LAT_all ) ) then
       ! unexpected
       ! INTERP_domain_compatibility should been called
       LOG_ERROR("get_IJrange",*) "unexpected error", LAT_min, LAT_max, minval( LAT_all ), maxval( LAT_all )
       call PRC_abort
    end if
    min = maxval( minval( LAT_all(:,:), dim=1 ), mask=all( LAT_all(:,:) < LAT_min, dim=1 ) )
    max = minval( maxval( LAT_all(:,:), dim=1 ), mask=all( LAT_all(:,:) > LAT_max, dim=1 ) )
    LAT_mask(:) = any( LAT_all(:,:) - min > -EPS, dim=1 ) .AND. any( LAT_all(:,:) - max < EPS, dim=1 )
    do j = 1, JA_org
       if( LAT_mask(j) ) then
          JS_org = j
          exit
       endif
    end do
    do j = JA_org, 1, -1
       if( LAT_mask(j) ) then
          JE_org = j
          exit
       endif
    end do

    return
  end subroutine get_IJrange

end module mod_realinput
