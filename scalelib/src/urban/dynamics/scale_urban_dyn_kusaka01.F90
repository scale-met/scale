!-------------------------------------------------------------------------------
!> module urban / dynamics / Kusaka01
!!
!! @par Description
!!          Single-layer Urban Canopy Model (Kusaka et al. 2001, BLM)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_urban_dyn_kusaka01
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_cpl_sfc_index
  use scale_mapprojection, only: &
      BASE_LON => MAPPROJECTION_basepoint_lon
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_DYN_kusaka01_setup
  public :: URBAN_DYN_kusaka01_finalize
  public :: URBAN_DYN_kusaka01

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: SLC_main
  private :: canopy_wind
  private :: cal_beta
  private :: cal_psi
  private :: mos
  private :: multi_layer
  private :: urban_param_setup
  private :: read_urban_param_table
  private :: read_urban_gridded_data_2D
  private :: read_urban_gridded_data_3D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  ! from namelist
  real(RP), private :: DTS_MAX    =    0.1_RP ! maximum dT during one step [K/step]
                                              ! DTS_MAX * dt
  integer,  private :: BOUND      =    1      ! Boundary Condition for Roof, Wall, Ground Layer Temp
                                              !       [1: Zero-Flux, 2: T = Constant]
#ifdef DEBUG_KUSAKA01
  logical,  private :: debug      = .true.
#else
  logical,  private :: debug      = .false.
#endif
  !$acc declare create(DTS_MAX,BOUND,debug)

  ! urban parameters
  real(RP), private :: ZR         =   10.0_RP ! roof level (building height) [m]
  real(RP), private :: roof_width =    9.0_RP ! roof width [m]
  real(RP), private :: road_width =   11.0_RP ! road width [m]
  real(RP), private :: SIGMA_ZED  =    1.0_RP ! Standard deviation of roof height [m]
  real(RP), private :: AH_TBL     =   17.5_RP ! Sensible Anthropogenic heat from urban subgrid [W/m^2]
  real(RP), private :: AHL_TBL    =    0.0_RP ! Latent Anthropogenic heat from urban subgrid [W/m^2]
  real(RP), private :: BETR_CONST =   -1.0_RP ! Evaporation efficiency of roof     [-]
  real(RP), private :: BETB_CONST =   -1.0_RP !                        of building [-]
  real(RP), private :: BETG_CONST =   -1.0_RP !                        of ground   [-]
  real(RP), private :: STRGR      =    0.0_RP ! rain strage on roof     [-]
  real(RP), private :: STRGB      =    0.0_RP !             on wall     [-]
  real(RP), private :: STRGG      =    0.0_RP !             on ground   [-]
  real(RP), private :: CAPR       =  1.2E6_RP ! heat capacity of roof   [J m-3 K]
  real(RP), private :: CAPB       =  1.2E6_RP !               of wall   [J m-3 K]
  real(RP), private :: CAPG       =  1.2E6_RP !               of ground [J m-3 K]
  real(RP), private :: AKSR       =   2.28_RP ! thermal conductivity of roof   [W m-1 K]
  real(RP), private :: AKSB       =   2.28_RP !                      of wall   [W m-1 K]
  real(RP), private :: AKSG       =   2.28_RP !                      of ground [W m-1 K]
  real(RP), private :: ALBR       =    0.2_RP ! surface albedo of roof
  real(RP), private :: ALBB       =    0.2_RP ! surface albedo of wall
  real(RP), private :: ALBG       =    0.2_RP ! surface albedo of ground
  real(RP), private :: EPSR       =   0.90_RP ! Surface emissivity of roof
  real(RP), private :: EPSB       =   0.90_RP ! Surface emissivity of wall
  real(RP), private :: EPSG       =   0.90_RP ! Surface emissivity of ground
  real(RP), private :: Z0R        =   0.01_RP ! roughness length for momentum of building roof
  real(RP), private :: Z0B        = 0.0001_RP ! roughness length for momentum of building wall
  real(RP), private :: Z0G        =   0.01_RP ! roughness length for momentum of ground
  real(RP), private :: TRLEND     = 293.00_RP ! lower boundary condition of roof temperature [K]
  real(RP), private :: TBLEND     = 293.00_RP ! lower boundary condition of wall temperature [K]
  real(RP), private :: TGLEND     = 293.00_RP ! lower boundary condition of ground temperature [K]
  !$acc declare create(ZR,BETR_CONST,BETB_CONST,BETG_CONST,STRGR,STRGB,STRGG,CAPR,CAPB,CAPG,AKSR,AKSB,AKSG,ALBR,ALBB,ALBG,EPSR,EPSB,EPSG,Z0R,TRLEND,TBLEND,TGLEND)

  ! calculated in subroutine urban_param_set
  real(RP), private :: R                       ! Normalized roof wight (eq. building coverage ratio)
  real(RP), private :: RW                      ! (= 1 - R)
  real(RP), private :: HGT                     ! Normalized building height
  real(RP), private :: Z0HR                    ! roughness length for heat of roof
  real(RP), private :: Z0HB                    ! roughness length for heat of building wall
  real(RP), private :: Z0HG                    ! roughness length for heat of ground
  real(RP), private :: Z0C_TBL                 ! Roughness length above canyon for momentum [m]
  real(RP), private :: Z0HC_TBL                ! Roughness length above canyon for heat [m]
  real(RP), private :: ZDC_TBL                 ! Displacement height [m]
  real(RP), private :: SVF                     ! Sky view factor [-]
  !$acc declare create(R,RW,HGT,Z0HR,Z0HB,Z0HG,SVF)

  ! history
  integer, private :: I_SHR, I_SHB, I_SHG
  integer, private :: I_LHR, I_LHB, I_LHG
  integer, private :: I_GHR, I_GHB, I_GHG
  integer, private :: I_RNR, I_RNB, I_RNG, I_RNgrd

#ifdef DEBUG_KUSAKA01
  integer, private :: cnt_num1, cnt_num2
  integer, private :: cnt_itr1, cnt_itr2
  integer, private :: max_itr1, max_itr2
  !$acc declare create(cnt_num1,cnt_num2,cnt_itr1,cnt_itr2,max_itr1,max_itr2)
#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_DYN_kusaka01_setup( &
       UIA, UIS, UIE, UJA, UJS, UJE,   &
       fact_urban,                     &
       Z0M, Z0H, Z0E, ZD,              &
       AH_URB, AHL_URB, AH_TOFFSET     )
    use scale_prc, only: &
       PRC_myrank,       &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_file_history, only: &
       FILE_HISTORY_reg
    implicit none

    integer,  intent(in)  :: UIA, UIS, UIE
    integer,  intent(in)  :: UJA, UJS, UJE
    real(RP), intent(in)  :: fact_urban(UIA,UJA)
    real(RP), intent(out) :: Z0M(UIA,UJA)
    real(RP), intent(out) :: Z0H(UIA,UJA)
    real(RP), intent(out) :: Z0E(UIA,UJA)
    real(RP), intent(out) :: ZD (UIA,UJA)
    real(RP), intent(out) :: AH_URB  (UIA,UJA,1:24)
    real(RP), intent(out) :: AHL_URB (UIA,UJA,1:24)
    real(RP), intent(out) :: AH_TOFFSET

    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME       = ''          !< urban parameter table
    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_Z0M_IN_FILENAME = ''          !< gridded data of Z0M
    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_Z0M_IN_VARNAME  = 'URBAN_Z0M' !< var name of gridded data for Z0M
    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_Z0H_IN_FILENAME = ''          !< gridded data of Z0H
    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_Z0H_IN_VARNAME  = 'URBAN_Z0H' !< var name of gridded data for Z0H
    !character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_ZD_IN_FILENAME  = ''          !< gridded data of ZD
    !character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_ZD_IN_VARNAME   = 'URBAN_ZD'  !< var name of gridded data for Zd
    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_AH_IN_FILENAME  = ''          !< gridded data of AH
    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_AH_IN_VARNAME   = 'URBAN_AH'  !< var name of gridded data for AH
    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_AHL_IN_FILENAME = ''          !< gridded data of AHL
    character(len=H_LONG) :: URBAN_DYN_KUSAKA01_GRIDDED_AHL_IN_VARNAME  = 'URBAN_AHL' !< var name of gridded data for AHL

    namelist / PARAM_URBAN_DYN_KUSAKA01 /          &
       DTS_MAX,                                    &
       BOUND,                                      &
       debug,                                      &
       URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME,       &
       URBAN_DYN_KUSAKA01_GRIDDED_Z0M_IN_FILENAME, &
       URBAN_DYN_KUSAKA01_GRIDDED_Z0H_IN_FILENAME, &
       !URBAN_DYN_KUSAKA01_GRIDDED_ZD_IN_FILENAME,  &
       URBAN_DYN_KUSAKA01_GRIDDED_AH_IN_FILENAME,  &
       URBAN_DYN_KUSAKA01_GRIDDED_AHL_IN_FILENAME

    real(RP) :: udata(UIA,UJA)
    real(RP) :: udata2(UIA,UJA,24)
    real(RP) :: rtime
    integer  :: itime

    real(RP) :: ahdiurnal(1:24)        ! AH diurnal profile

    integer  :: i, j, k
    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_DYN_kusaka01_setup",*) 'Setup'

    !$acc data copyin(fact_urban) &
    !$acc      copyout(Z0M, Z0H, Z0E, ZD, AH_URB, AHL_URB)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_DYN_KUSAKA01,iostat=ierr)
    if( ierr < 0 ) then     !--- missing
       LOG_INFO("URBAN_DYN_kusaka01_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("URBAN_DYN_kusaka01_setup",*) 'Not appropriate names in namelist PARAM_URBAN_DYN_KUSAKA01. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_URBAN_DYN_KUSAKA01)

    !-- read urban parameter from file
    if( URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME /= '' ) then
     call read_urban_param_table( trim(URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME) )
    endif

    ! set other urban parameters
    call urban_param_setup

    ! Local time: 1 to 24
    ahdiurnal(:) = (/ 0.356, 0.274, 0.232, 0.251, 0.375, 0.647, 0.919, 1.135, 1.249, 1.328, &
                      1.365, 1.363, 1.375, 1.404, 1.457, 1.526, 1.557, 1.521, 1.372, 1.206, &
                      1.017, 0.876, 0.684, 0.512                                            /)

    ! Shift to UTC based on local solar timezone of domain center
    rtime = modulo(BASE_LON, 360.0_RP) / 15.0_RP
    itime = nint(rtime)
    AH_TOFFSET = rtime - itime                 ! currently not used: difference of time from domain center
    ahdiurnal(:) = cshift(ahdiurnal(:), -1*itime) ! convert from LT to UTC

    do j = UJS, UJE
    do i = UIS, UIE
       Z0M(i,j) = Z0C_TBL
       Z0H(i,j) = Z0HC_TBL
       Z0E(i,j) = Z0HC_TBL
       ZD(i,j)  = ZDC_TBL
       do k = 1, 24 ! UTC
          AH_URB (i,j,k) = AH_TBL  * ahdiurnal(k) * fact_urban(i,j)
          AHL_URB(i,j,k) = AHL_TBL * ahdiurnal(k) * fact_urban(i,j)
       enddo
    enddo
    enddo

    !-- replace gridded Z0M data if there is a file
    if( URBAN_DYN_KUSAKA01_GRIDDED_Z0M_IN_FILENAME /= '' ) then
     udata(:,:) = 0.0_RP
     call read_urban_gridded_data_2D(                   &
            UIA, UJA,                                   &
            URBAN_DYN_KUSAKA01_GRIDDED_Z0M_IN_FILENAME, &
            URBAN_DYN_KUSAKA01_GRIDDED_Z0M_IN_VARNAME,  &
            udata(:,:)                                  )

      ! replace to gridded data
      do j = UJS, UJE
      do i = UIS, UIE
         if( udata(i,j) /= UNDEF )then
            if ( udata(i,j) > 0.0_RP ) then
               Z0M(i,j) = udata(i,j)
            else if ( udata(i,j) < 0.0_RP ) then
               LOG_ERROR("URBAN_DYN_kusaka01_setup",*) 'Gridded Z0M data includes data less than 0. Please check data!',i,j
               call PRC_abort
            else ! Z0M = 0[m]
               LOG_WARN("URBAN_DYN_kusaka01_setup",*) 'Gridded Z0M data includes 0; default or table value is used to avoid zero division',PRC_myrank,i,j
            endif
         endif
      enddo
      enddo
    endif

    !-- read gridded Z0H & Z0E data from a file
    if( URBAN_DYN_KUSAKA01_GRIDDED_Z0H_IN_FILENAME /= '' ) then
     udata = 0.0_RP
     call read_urban_gridded_data_2D(                   &
            UIA, UJA,                                   &
            URBAN_DYN_KUSAKA01_GRIDDED_Z0H_IN_FILENAME, &
            URBAN_DYN_KUSAKA01_GRIDDED_Z0H_IN_VARNAME,  &
            udata                                       )

      ! replace to gridded data
      do j = UJS, UJE
      do i = UIS, UIE
         if( udata(i,j) /= UNDEF )then
            if ( udata(i,j) > 0.0_RP ) then
               Z0H(i,j) = udata(i,j)
               Z0E(i,j) = udata(i,j)
            else if ( udata(i,j) < 0.0_RP ) then
               LOG_ERROR("URBAN_DYN_kusaka01_setup",*) 'Gridded Z0H data includes data less than 0. Please check data!',i,j
               call PRC_abort
            else ! Z0H = 0[m]
               LOG_WARN("URBAN_DYN_kusaka01_setup",*) 'Gridded Z0H data includes 0; default or table value is used to avoid zero division',PRC_myrank,i,j
            endif
         endif
      enddo
      enddo
    endif

    ! currently NOT USED
    !-- read gridded ZD data from a file
    !if( URBAN_DYN_KUSAKA01_GRIDDED_ZD_IN_FILENAME /= '' ) then
    ! udata = 0.0_RP
    ! call read_urban_gridded_data_2D(                  &
    !        UIA, UJA,                                  &
    !        URBAN_DYN_KUSAKA01_GRIDDED_ZD_IN_FILENAME, &
    !        URBAN_DYN_KUSAKA01_GRIDDED_ZD_IN_VARNAME,  &
    !        udata                                      )
    !
    !  ! replace to gridded data
    !  do j = UJS, UJE
    !  do i = UIS, UIE
    !     if( udata(i,j) /= UNDEF )then
    !        if ( udata(i,j) >= 0.0_RP ) then
    !           ZD(i,j) = udata(i,j)
    !        else
    !           LOG_ERROR("URBAN_DYN_kusaka01_setup",*) 'Gridded ZD data includes data less than 0. Please check data!',i,j
    !           call PRC_abort
    !        endif
    !     endif
    !  enddo
    !  enddo
    !endif

    !-- read gridded AH data from a file
    if( URBAN_DYN_KUSAKA01_GRIDDED_AH_IN_FILENAME /= '' ) then
       udata2 = 0.0_RP
       call read_urban_gridded_data_3D(                  &
            UIA, UJA,                                  &
            URBAN_DYN_KUSAKA01_GRIDDED_AH_IN_FILENAME, &
            URBAN_DYN_KUSAKA01_GRIDDED_AH_IN_VARNAME,  &
            udata2                                     )

       ! replace to gridded data
       do k = 1, 24
       do j = UJS, UJE
       do i = UIS, UIE
          if( udata2(i,j,k) /= UNDEF )then
             AH_URB(i,j,k) = udata2(i,j,k)
          endif
       enddo
       enddo
       enddo
    endif

    !-- read gridded AHL data from a file
    if( URBAN_DYN_KUSAKA01_GRIDDED_AHL_IN_FILENAME /= '' ) then
     udata2 = 0.0_RP
     call read_urban_gridded_data_3D(                   &
            UIA, UJA,                                   &
            URBAN_DYN_KUSAKA01_GRIDDED_AHL_IN_FILENAME, &
            URBAN_DYN_KUSAKA01_GRIDDED_AHL_IN_VARNAME,  &
            udata2                                      )

      ! replace to gridded data
      do k = 1, 24
      do j = UJS, UJE
      do i = UIS, UIE
         if( udata2(i,j,k) /= UNDEF )then
            AHL_URB(i,j,k) = udata2(i,j,k)
         endif
      enddo
      enddo
      enddo
    endif

    call FILE_HISTORY_reg( 'URBAN_SHR',   'urban sensible heat flux on roof',    'W/m2', I_SHR  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_SHB',   'urban sensible heat flux on wall',    'W/m2', I_SHB  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_SHG',   'urban sensible heat flux on road',    'W/m2', I_SHG  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_LHR',   'urban latent heat flux on roof',      'W/m2', I_LHR  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_LHB',   'urban latent heat flux on wall',      'W/m2', I_LHB  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_LHG',   'urban latent heat flux on road',      'W/m2', I_LHG  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_GHR',   'urban ground heat flux on roof',      'W/m2', I_GHR  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_GHB',   'urban ground heat flux on wall',      'W/m2', I_GHB  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_GHG',   'urban ground heat flux on road',      'W/m2', I_GHG  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_RNR',   'urban net radiation on roof',         'W/m2', I_RNR  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_RNB',   'urban net radiation on wall',         'W/m2', I_RNB  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_RNG',   'urban net radiation on road',         'W/m2', I_RNG  , ndims=2 )
    call FILE_HISTORY_reg( 'URBAN_RNgrd', 'urban grid average of net radiation', 'W/m2', I_RNgrd, ndims=2 )

    !$acc update device(Z0M, Z0H, Z0E, ZD, AH_URB, AHL_URB)
    !$acc end data

    !$acc update device(DTS_MAX,BOUND,debug)
    !$acc update device(ZR,BETR_CONST,BETB_CONST,BETG_CONST,STRGR,STRGB,STRGG,CAPR,CAPB,CAPG,AKSR,AKSB,AKSG,ALBR,ALBB,ALBG,EPSR,EPSB,EPSG,Z0R,TRLEND,TBLEND,TGLEND)
    !$acc update device(R,RW,HGT,Z0HR,Z0HB,Z0HG,SVF)

#ifdef DEBUG_KUSAKA01
    cnt_num1 = 0
    cnt_num2 = 0
    cnt_itr1 = 0
    cnt_itr2 = 0
    max_itr1 = 0
    max_itr2 = 0
    !$acc update device(cnt_num1,cnt_num2,cnt_itr1,cnt_itr2,max_itr1,max_itr2)
#endif

    return
  end subroutine URBAN_DYN_kusaka01_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine URBAN_DYN_kusaka01_finalize
    use mpi
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_IsMaster
    implicit none

#ifdef DEBUG_KUSAKA01
    integer :: iwork1(4), iwork2(2), ierr
    !$acc update host(cnt_num1,cnt_num2,cnt_itr1,cnt_itr2,max_itr1,max_itr2)
    iwork1(1) = cnt_num1
    iwork1(2) = cnt_num2
    iwork1(3) = cnt_itr1
    iwork1(4) = cnt_itr2
    call MPI_Reduce( MPI_IN_PLACE, iwork1, 4, MPI_INTEGER, MPI_SUM, 0, PRC_LOCAL_COMM_WORLD, ierr)
    iwork2(1) = max_itr1
    iwork2(2) = max_itr2
    call MPI_Reduce( MPI_IN_PLACE, iwork2, 2, MPI_INTEGER, MPI_MAX, 0, PRC_LOCAL_COMM_WORLD, ierr)

    if ( PRC_IsMaster ) then
       LOG_INFO("URBAN_DYN_kusaka01_finalize",*) "Averaged iteration count"
       LOG_INFO_CONT(*) "TR:    ", real(iwork1(3),DP) / iwork1(1), ", (max ", iwork2(1), ")"
       LOG_INFO_CONT(*) "TB/TG: ", real(iwork1(4),DP) / iwork1(2), ", (max ", iwork2(2), ")"
    end if
#endif

  end subroutine URBAN_DYN_kusaka01_finalize

  !-----------------------------------------------------------------------------
  !> Main routine for land submodel

  subroutine URBAN_DYN_kusaka01( &
       UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
       TMPA, PRSA,                      &
       U1, V1,                          &
       DENS, QA, LHV,                   &
       Z1,                              &
       RHOS, PRSS,                      &
       LWD, SWD,                        &
       RAIN, EFLX,                      &
       Z0M, Z0H, Z0E,                   &
       ZD,                              &
       CDZ,                             &
       TanSL_X, TanSL_Y,                &
       fact_urban,                      &
       dt,                              &
       TRL_URB, TBL_URB, TGL_URB,       &
       TR_URB, TB_URB, TG_URB,          &
       TC_URB, QC_URB, UC_URB,          &
       RAINR_URB, RAINB_URB, RAING_URB, &
       ROFF_URB,                        &
       SFC_TEMP,                        &
       ALBEDO,                          &
       MWFLX, MUFLX, MVFLX,             &
       SHFLX, LHFLX, GHFLX,             &
       Ustar, Tstar, Qstar, Wstar,      &
       RLmo,                            &
       U10, V10, T2, Q2                 )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       EPS   => CONST_EPS,  &
       UNDEF => CONST_UNDEF, &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap
    use scale_bulkflux, only: &
       BULKFLUX
    implicit none
    integer, intent(in) :: UKA, UKS, UKE
    integer, intent(in) :: UIA, UIS, UIE
    integer, intent(in) :: UJA, UJS, UJE

    real(RP), intent(in) :: TMPA(UIA,UJA)
    real(RP), intent(in) :: PRSA(UIA,UJA)
    real(RP), intent(in) :: U1  (UIA,UJA)
    real(RP), intent(in) :: V1  (UIA,UJA)
    real(RP), intent(in) :: DENS(UIA,UJA)
    real(RP), intent(in) :: QA  (UIA,UJA)
    real(RP), intent(in) :: LHV (UIA,UJA)
    real(RP), intent(in) :: Z1  (UIA,UJA)
    real(RP), intent(in) :: RHOS(UIA,UJA) ! density  at the surface [kg/m3]
    real(RP), intent(in) :: PRSS(UIA,UJA)
    real(RP), intent(in) :: LWD (UIA,UJA,2)
    real(RP), intent(in) :: SWD (UIA,UJA,2)
    real(RP), intent(in) :: RAIN(UIA,UJA)
    real(RP), intent(in) :: EFLX(UIA,UJA)
    real(RP), intent(in) :: Z0M      (UIA,UJA)
    real(RP), intent(in) :: Z0H      (UIA,UJA)
    real(RP), intent(in) :: Z0E      (UIA,UJA)
    real(RP), intent(in) :: ZD       (UIA,UJA)
    real(RP), intent(in) :: CDZ(UKA)
    real(RP), intent(in) :: TanSL_X(UIA,UJA)
    real(RP), intent(in) :: TanSL_Y(UIA,UJA)
    real(RP), intent(in) :: fact_urban(UIA,UJA)
    real(DP), intent(in) :: dt

    real(RP), intent(inout) :: TR_URB   (UIA,UJA)
    real(RP), intent(inout) :: TB_URB   (UIA,UJA)
    real(RP), intent(inout) :: TG_URB   (UIA,UJA)
    real(RP), intent(inout) :: TC_URB   (UIA,UJA)
    real(RP), intent(inout) :: QC_URB   (UIA,UJA)
    real(RP), intent(inout) :: UC_URB   (UIA,UJA)
    real(RP), intent(inout) :: TRL_URB  (UKS:UKE,UIA,UJA)
    real(RP), intent(inout) :: TBL_URB  (UKS:UKE,UIA,UJA)
    real(RP), intent(inout) :: TGL_URB  (UKS:UKE,UIA,UJA)
    real(RP), intent(inout) :: RAINR_URB(UIA,UJA)
    real(RP), intent(inout) :: RAINB_URB(UIA,UJA)
    real(RP), intent(inout) :: RAING_URB(UIA,UJA)
    real(RP), intent(out)   :: ROFF_URB (UIA,UJA)

    real(RP), intent(out) :: SFC_TEMP(UIA,UJA)
    real(RP), intent(out) :: ALBEDO  (UIA,UJA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(out) :: MWFLX   (UIA,UJA)
    real(RP), intent(out) :: MUFLX   (UIA,UJA)
    real(RP), intent(out) :: MVFLX   (UIA,UJA)
    real(RP), intent(out) :: SHFLX   (UIA,UJA)
    real(RP), intent(out) :: LHFLX   (UIA,UJA)
    real(RP), intent(out) :: GHFLX   (UIA,UJA)
    real(RP), intent(out) :: Ustar   (UIA,UJA)
    real(RP), intent(out) :: Tstar   (UIA,UJA)
    real(RP), intent(out) :: Qstar   (UIA,UJA)
    real(RP), intent(out) :: Wstar   (UIA,UJA)
    real(RP), intent(out) :: RLmo    (UIA,UJA)
    real(RP), intent(out) :: U10     (UIA,UJA)
    real(RP), intent(out) :: V10     (UIA,UJA)
    real(RP), intent(out) :: T2      (UIA,UJA)
    real(RP), intent(out) :: Q2      (UIA,UJA)

    ! parameter
    logical,  parameter :: LSOLAR = .false. ! [true=both, false=SSG only]

    real(RP), parameter :: Uabs_min = 0.1_RP

    ! work
    real(RP) :: TR
    real(RP) :: TB
    real(RP) :: TG
    real(RP) :: TC
    real(RP) :: QC
    real(RP) :: UC
    real(RP) :: TRL(UKS:UKE)
    real(RP) :: TBL(UKS:UKE)
    real(RP) :: TGL(UKS:UKE)
    real(RP) :: RAINR
    real(RP) :: RAINB
    real(RP) :: RAING
    real(RP) :: ALBD_LW
    real(RP) :: ALBD_SW

    real(RP) :: SHR  (UIA,UJA)
    real(RP) :: SHB  (UIA,UJA)
    real(RP) :: SHG  (UIA,UJA)
    real(RP) :: LHR  (UIA,UJA)
    real(RP) :: LHB  (UIA,UJA)
    real(RP) :: LHG  (UIA,UJA)
    real(RP) :: GHR  (UIA,UJA)
    real(RP) :: GHB  (UIA,UJA)
    real(RP) :: GHG  (UIA,UJA)
    real(RP) :: RNR  (UIA,UJA)
    real(RP) :: RNB  (UIA,UJA)
    real(RP) :: RNG  (UIA,UJA)
    real(RP) :: RNgrd(UIA,UJA)

    real(RP) :: DZR(UKA)     ! thickness of each roof layer [m]
    real(RP) :: DZB(UKA)     ! thickness of each building layer [m]
    real(RP) :: DZG(UKA)     ! thickness of each road layer [m]

    real(RP) :: SWDt(2)
    real(RP) :: LWDt(2)

    real(RP) :: Uabs  ! modified absolute velocity [m/s]

    real(RP) :: MFLX
    real(RP) :: w

    ! work
    real(RP) :: TRLP(UKA)
    real(RP) :: TBLP(UKA)
    real(RP) :: TGLP(UKA)
    real(RP) :: A(UKA)
    real(RP) :: B(UKA)
    real(RP) :: C(UKA)
    real(RP) :: D(UKA)
    real(RP) :: P(UKA)
    real(RP) :: Q(UKA)

    logical :: converged

#ifdef DEBUG_KUSAKA01
    integer :: cn1, cn2, ci1, ci2, mi1, mi2
#endif

    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'urban / dynamics / Kusaka01'

    !$acc data copyin(TMPA,PRSA,U1,V1,DENS,QA,LHV,Z1,RHOS,PRSS,LWD,SWD,RAIN,EFLX,Z0M,Z0H,Z0E,ZD,CDZ,TanSL_X,TanSL_Y,fact_urban) &
    !$acc      copy(TR_URB,TB_URB,TG_URB,TC_URB,QC_URB,UC_URB,TRL_URB,TBL_URB,TGL_URB,RAINR_URB,RAINB_URB,RAING_URB) &
    !$acc      copyout(ROFF_URB,SFC_TEMP,ALBEDO,MWFLX,MUFLX,MVFLX,SHFLX,LHFLX,GHFLX,Ustar,Tstar,Qstar,Wstar,RLmo,U10,V10,T2,Q2) &
    !$acc      create(SHR,SHB,SHG,LHR,LHB,LHG,GHR,GHB,GHG,RNR,RNB,RNG,RNgrd,DZR,DZB,DZG)


    !$acc kernels
    DZR(:) = CDZ(:)
    DZB(:) = CDZ(:)
    DZG(:) = CDZ(:)
    !$acc end kernels

    converged = .true.

#ifdef DEBUG_KUSAKA01
    cn1 = 0; cn2 = 0; ci1 = 0; ci2 = 0; mi1 = 0; mi2 = 0
#endif

    !$omp parallel do schedule(dynamic) collapse(2) &
#ifdef DEBUG_KUSAKA01
    !$omp reduction(+:cn1,cn2,ci1,ci2) reduction(max:mi1,mi2) &
#endif
    !$omp private(w,Uabs,TR,TB,TG,TC,QC,UC,TRL,TBL,TGL,RAINR,RAINB,RAING,ALBD_LW,ALBD_SW,MFLX,SWDt,LWDt, &
    !$omp         TRLP,TBLP,TGLP,A,B,C,D,P,Q)
    !$acc parallel
    !$acc loop gang collapse(2) reduction(.and.: converged) independent &
#ifdef DEBUG_KUSAKA01
    !$acc reduction(+:cn1,cn2,ci1,ci2) reduction(max:mi1,mi2) &
#endif
    !$acc private(w,Uabs,TR,TB,TG,TC,QC,UC,TRL,TBL,TGL,RAINR,RAINB,RAING,ALBD_LW,ALBD_SW,MFLX,SWDt,LWDt, &
    !$acc         TRLP,TBLP,TGLP,A,B,C,D,P,Q)
    do j = UJS, UJE
    do i = UIS, UIE

    if( fact_urban(i,j) > 0.0_RP ) then

       w = U1(i,j) * TanSL_X(i,j) + V1(i,j) * TanSL_Y(i,j)
       Uabs = sqrt( U1(i,j)**2 + V1(i,j)**2 + w**2 )

       ! save
       TR = TR_URB(i,j)
       TB = TB_URB(i,j)
       TG = TG_URB(i,j)
       TC = TC_URB(i,j)
       QC = QC_URB(i,j)
       UC = UC_URB(i,j)

       do k = UKS, UKE
          TRL(k) = TRL_URB(k,i,j)
          TBL(k) = TBL_URB(k,i,j)
          TGL(k) = TGL_URB(k,i,j)
       end do

       RAINR = RAINR_URB(i,j)
       RAINB = RAINB_URB(i,j)
       RAING = RAING_URB(i,j)

       SWDt(:) = SWD(i,j,:)
       LWDt(:) = LWD(i,j,:)

       call SLC_main( UKA, UKS, UKE, &
                      TRL     (:),        & ! [INOUT]
                      TBL     (:),        & ! [INOUT]
                      TGL     (:),        & ! [INOUT]
                      TR,                 & ! [INOUT]
                      TB,                 & ! [INOUT]
                      TG,                 & ! [INOUT]
                      TC,                 & ! [INOUT]
                      QC,                 & ! [INOUT]
                      UC,                 & ! [INOUT]
                      RAINR,              & ! [INOUT]
                      RAINB,              & ! [INOUT]
                      RAING,              & ! [INOUT]
                      ROFF_URB(i,j),      & ! [OUT]
                      ALBD_LW,            & ! [OUT]
                      ALBD_SW,            & ! [OUT]
                      SHR     (i,j),      & ! [OUT]
                      SHB     (i,j),      & ! [OUT]
                      SHG     (i,j),      & ! [OUT]
                      LHR     (i,j),      & ! [OUT]
                      LHB     (i,j),      & ! [OUT]
                      LHG     (i,j),      & ! [OUT]
                      GHR     (i,j),      & ! [OUT]
                      GHB     (i,j),      & ! [OUT]
                      GHG     (i,j),      & ! [OUT]
                      RNR     (i,j),      & ! [OUT]
                      RNB     (i,j),      & ! [OUT]
                      RNG     (i,j),      & ! [OUT]
                      SFC_TEMP(i,j),      & ! [OUT]
                      RNgrd   (i,j),      & ! [OUT]
                      SHFLX   (i,j),      & ! [OUT]
                      LHFLX   (i,j),      & ! [OUT]
                      GHFLX   (i,j),      & ! [OUT]
                      MFLX,               & ! [OUT]
                      Ustar   (i,j),      & ! [OUT]
                      Tstar   (i,j),      & ! [OUT]
                      Qstar   (i,j),      & ! [OUT]
                      U10     (i,j),      & ! [OUT]
                      V10     (i,j),      & ! [OUT]
                      T2      (i,j),      & ! [OUT]
                      Q2      (i,j),      & ! [OUT]
                      converged,          & ! [OUT]
                      LSOLAR,             & ! [IN]
                      PRSA    (i,j),      & ! [IN]
                      PRSS    (i,j),      & ! [IN]
                      TMPA    (i,j),      & ! [IN]
                      RHOS    (i,j),      & ! [IN]
                      QA      (i,j),      & ! [IN]
                      Uabs,               & ! [IN]
                      U1      (i,j),      & ! [IN]
                      V1      (i,j),      & ! [IN]
                      LHV     (i,j),      & ! [IN]
                      Z1      (i,j),      & ! [IN]
                      SWDt    (:),        & ! [IN]
                      LWDt    (:),        & ! [IN]
                      RAIN    (i,j),      & ! [IN]
                      EFLX    (i,j),      & ! [IN]
                      DENS    (i,j),      & ! [IN]
                      Z0M     (i,j),      & ! [IN]
                      Z0H     (i,j),      & ! [IN]
                      ZD      (i,j),      & ! [IN]
                      DZR(:),             & ! [IN]
                      DZG(:), DZB(:),     & ! [IN]
                      dt,                 & ! [IN]
                      i, j,               & ! [IN]
#ifdef DEBUG_KUSAKA01
                      cn1, cn2, ci1, ci2, mi1, mi2, & ! [INOUT]
#endif
                      TRLP(:), TBLP(:), TGLP(:), & ! (work)
                      A(:), B(:), C(:), D(:),    & ! (work)
                      P(:), Q(:)                 ) ! (work)

       ! update
       TR_URB(i,j) = TR
       TB_URB(i,j) = TB
       TG_URB(i,j) = TG
       TC_URB(i,j) = TC
       QC_URB(i,j) = QC
       UC_URB(i,j) = UC
       do k = UKS, UKE
          TRL_URB(k,i,j) = TRL(k)
          TBL_URB(k,i,j) = TBL(k)
          TGL_URB(k,i,j) = TGL(k)
       end do
       RAINR_URB(i,j) = RAINR
       RAINB_URB(i,j) = RAINB
       RAING_URB(i,j) = RAING

       ALBEDO(i,j,I_R_direct ,I_R_IR ) = ALBD_LW
       ALBEDO(i,j,I_R_diffuse,I_R_IR ) = ALBD_LW
       ALBEDO(i,j,I_R_direct ,I_R_NIR) = ALBD_SW
       ALBEDO(i,j,I_R_diffuse,I_R_NIR) = ALBD_SW
       ALBEDO(i,j,I_R_direct ,I_R_VIS) = ALBD_SW
       ALBEDO(i,j,I_R_diffuse,I_R_VIS) = ALBD_SW

       if ( Uabs < EPS ) then
          MWFLX(i,j) = 0.0_RP
          MUFLX(i,j) = 0.0_RP
          MVFLX(i,j) = 0.0_RP
       else
          MWFLX(i,j) = MFLX * w / Uabs
          MUFLX(i,j) = MFLX * U1(i,j) / Uabs
          MVFLX(i,j) = MFLX * V1(i,j) / Uabs
       end if

    else
       SFC_TEMP(i,j)     = UNDEF
       ALBEDO  (i,j,:,:) = UNDEF
       MWFLX   (i,j)     = UNDEF
       MUFLX   (i,j)     = UNDEF
       MVFLX   (i,j)     = UNDEF
       SHFLX   (i,j)     = UNDEF
       LHFLX   (i,j)     = UNDEF
       GHFLX   (i,j)     = UNDEF
       Ustar   (i,j)     = UNDEF
       Tstar   (i,j)     = UNDEF
       Qstar   (i,j)     = UNDEF
       Wstar   (i,j)     = UNDEF
       RLmo    (i,j)     = UNDEF
       U10     (i,j)     = UNDEF
       V10     (i,j)     = UNDEF
       T2      (i,j)     = UNDEF
       Q2      (i,j)     = UNDEF
       SHR     (i,j)     = UNDEF
       SHB     (i,j)     = UNDEF
       SHG     (i,j)     = UNDEF
       LHR     (i,j)     = UNDEF
       LHB     (i,j)     = UNDEF
       LHG     (i,j)     = UNDEF
       GHR     (i,j)     = UNDEF
       GHB     (i,j)     = UNDEF
       GHG     (i,j)     = UNDEF
       RNR     (i,j)     = UNDEF
       RNB     (i,j)     = UNDEF
       RNG     (i,j)     = UNDEF
       RNgrd   (i,j)     = UNDEF
    endif

    end do
    end do
    !$acc end parallel

    if ( .not. converged ) then
       LOG_ERROR("URBAN_DYN_kusaka01_SLC_main",*) "not converged"
       call PRC_abort
    end if

    call put_history( UIA, UJA, &
                      SHR(:,:), SHB(:,:), SHG(:,:), &
                      LHR(:,:), LHB(:,:), LHG(:,:), &
                      GHR(:,:), GHB(:,:), GHG(:,:), &
                      RNR(:,:), RNB(:,:), RNG(:,:), &
                      RNgrd(:,:)                    )

    !$acc end data

#ifdef DEBUG_KUSAKA01
    cnt_num1 = cnt_num1 + cn1
    cnt_num2 = cnt_num2 + cn2
    cnt_itr1 = cnt_itr1 + ci1
    cnt_itr2 = cnt_itr2 + ci2
    max_itr1 = max(max_itr1, mi1)
    max_itr2 = max(max_itr2, mi2)
#endif

    return
  end subroutine URBAN_DYN_kusaka01

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine SLC_main( &
        UKA, UKS, UKE, &
        TRL,          & ! (inout)
        TBL,          & ! (inout)
        TGL,          & ! (inout)
        TR,           & ! (inout)
        TB,           & ! (inout)
        TG,           & ! (inout)
        TC,           & ! (inout)
        QC,           & ! (inout)
        UC,           & ! (inout)
        RAINR,        & ! (inout)
        RAINB,        & ! (inout)
        RAING,        & ! (inout)
        ROFF,         & ! (out)
        ALBD_LW_grid, & ! (out)
        ALBD_SW_grid, & ! (out)
        SHR,          & ! (out)
        SHB,          & ! (out)
        SHG,          & ! (out)
        LHR,          & ! (out)
        LHB,          & ! (out)
        LHG,          & ! (out)
        GHR,          & ! (out)
        GHB,          & ! (out)
        GHG,          & ! (out)
        RNR,          & ! (out)
        RNB,          & ! (out)
        RNG,          & ! (out)
        RTS,          & ! (out)
        RN,           & ! (out)
        SH,           & ! (out)
        LH,           & ! (out)
        GHFLX,        & ! (out)
        MFLX,         & ! (out)
        Ustar,        & ! (out)
        Tstar,        & ! (out)
        Qstar,        & ! (out)
        U10,          & ! (out)
        V10,          & ! (out)
        T2,           & ! (out)
        Q2,           & ! (out)
        converged,    & ! (out)
        LSOLAR,       & ! (in)
        PRSA,         & ! (in)
        PRSS,         & ! (in)
        TA,           & ! (in)
        RHOS,         & ! (in)
        QA,           & ! (in)
        UA,           & ! (in)
        U1,           & ! (in)
        V1,           & ! (in)
        LHV,          & ! (in)
        ZA,           & ! (in)
        SSG,          & ! (in)
        LLG,          & ! (in)
        RAIN,         & ! (in)
        EFLX,         & ! (in)
        RHOO,         & ! (in)
        Z0C,          & ! (in)
        Z0HC,         & ! (in)
        ZDC,          & ! (in)
        DZR, DZB, DZG, & ! (in)
        dt,           & ! (in)
        i, j,         & ! (in)
#ifdef DEBUG_KUSAKA01
        cnt_num1, cnt_num2, & ! (inout)
        cnt_itr1, cnt_itr2, & ! (inout)
        max_itr1, max_itr2, & ! (inout)
#endif
        TRLP, TBLP, TGLP, &
        A, B, C, D, P, Q )
    !$acc routine seq
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       EPS    => CONST_EPS,     &    ! small number (machine epsilon)
       KARMAN => CONST_KARMAN,  &    ! Kalman constant  [-]
       CPdry  => CONST_CPdry,   &    ! Heat capacity of dry air [J/K/kg]
       GRAV   => CONST_GRAV,    &    ! gravitational constant [m/s2]
       Rdry   => CONST_Rdry,    &    ! specific gas constant (dry) [J/kg/K]
       Rvap   => CONST_Rvap,    &    ! gas constant (water vapor) [J/kg/K]
       STB    => CONST_STB,     &    ! stefan-Boltzman constant [MKS unit]
       TEM00  => CONST_TEM00,   &    ! temperature reference (0 degree C) [K]
       PRE00  => CONST_PRE00         ! pressure reference [Pa]
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_atmos_saturation, only: &
       qsat => ATMOS_SATURATION_dens2qsat_liq, &
       dqs_dtem => ATMOS_SATURATION_dqs_dtem_dens_liq
    use scale_bulkflux, only: &
       BULKFLUX_diagnose_surface
    implicit none

    integer, intent(in) :: UKA, UKS, UKE

    !-- In/Out variables from/to Coupler to/from Urban
    real(RP), intent(inout) :: TRL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: TBL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: TGL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: TR   ! roof temperature              [K]
    real(RP), intent(inout) :: TB   ! building wall temperature     [K]
    real(RP), intent(inout) :: TG   ! road temperature              [K]
    real(RP), intent(inout) :: TC   ! urban-canopy air temperature  [K]
    real(RP), intent(inout) :: QC   ! urban-canopy air specific humidity [kg/kg]
    real(RP), intent(inout) :: UC   ! diagnostic canopy wind        [m/s]
    real(RP), intent(inout) :: RAINR ! rain amount in storage on roof     [kg/m2]
    real(RP), intent(inout) :: RAINB ! rain amount in storage on building [kg/m2]
    real(RP), intent(inout) :: RAING ! rain amount in storage on road     [kg/m2]
    real(RP), intent(out)   :: ROFF  ! runoff from urban           [kg/m2/s]

    !-- Output variables from Urban to Coupler
    real(RP), intent(out)   :: ALBD_SW_grid  ! grid mean of surface albedo for SW
    real(RP), intent(out)   :: ALBD_LW_grid  ! grid mean of surface albedo for LW ( 1-emiss )
    real(RP), intent(out)   :: SHR, SHB, SHG
    real(RP), intent(out)   :: LHR, LHB, LHG
    real(RP), intent(out)   :: GHR, GHB, GHG
    real(RP), intent(out)   :: RNR, RNB, RNG
    real(RP), intent(out)   :: RTS    ! radiative surface temperature    [K]
    real(RP), intent(out)   :: RN     ! net radition                     [W/m/m]
    real(RP), intent(out)   :: SH     ! sensible heat flux               [W/m/m]
    real(RP), intent(out)   :: LH     ! latent heat flux                 [W/m/m]
    real(RP), intent(out)   :: GHFLX  ! heat flux into the ground        [W/m/m]
    real(RP), intent(out)   :: MFLX   ! momentum flux                    [kg/m/s2]
    real(RP), intent(out)   :: Ustar  ! friction velocity                [m/s]
    real(RP), intent(out)   :: Tstar  ! temperature scale                [K]
    real(RP), intent(out)   :: Qstar  ! humidity scale                   [kg/kg]
    real(RP), intent(out)   :: U10    ! U wind at 10m                    [m/s]
    real(RP), intent(out)   :: V10    ! V wind at 10m                    [m/s]
    real(RP), intent(out)   :: T2     ! air temperature at 2m            [K]
    real(RP), intent(out)   :: Q2     ! specific humidity at 2m          [kg/kg]
    logical,  intent(out)   :: converged

    !-- configuration variables
    logical , intent(in)    :: LSOLAR ! logical   [true=both, false=SSG only]

    !-- Input variables from Coupler to Urban
    real(RP), intent(in)    :: PRSA ! Pressure at 1st atmospheric layer      [Pa]
    real(RP), intent(in)    :: PRSS ! Surface Pressure                       [Pa]
    real(RP), intent(in)    :: TA   ! temp at 1st atmospheric level          [K]
    real(RP), intent(in)    :: RHOS ! surface density                        [kg/m^3]
    real(RP), intent(in)    :: QA   ! specific humidity at 1st atmospheric level  [kg/kg]
    real(RP), intent(in)    :: UA   ! wind speed at 1st atmospheric level    [m/s]
    real(RP), intent(in)    :: U1   ! u at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: V1   ! v at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: LHV  ! latent heat of vaporization [J/kg]
    real(RP), intent(in)    :: ZA   ! height of 1st atmospheric level        [m]
    real(RP), intent(in)    :: SSG(2) ! downward total short wave radiation  [W/m/m]
    real(RP), intent(in)    :: LLG(2) ! downward long wave radiation         [W/m/m]
    real(RP), intent(in)    :: RAIN ! water flux                             [kg/m2/s]
    real(RP), intent(in)    :: EFLX ! internal energy flux                   [J/m2/s]
    real(RP), intent(in)    :: RHOO ! air density                            [kg/m^3]
    real(RP), intent(in)    :: Z0C  ! Roughness length above canyon for momentum [m]
    real(RP), intent(in)    :: Z0HC ! Roughness length above canyon for heat [m]
    real(RP), intent(in)    :: ZDC  ! Displacement height                    [m]
    real(RP), intent(in)    :: DZR(UKA)
    real(RP), intent(in)    :: DZB(UKA)
    real(RP), intent(in)    :: DZG(UKA)
    real(DP), intent(in)    :: dt

    integer,  intent(in)    :: i, j

#ifdef DEBUG_KUSAKA01
    integer, intent(inout) :: cnt_num1, cnt_num2
    integer, intent(inout) :: cnt_itr1, cnt_itr2
    integer, intent(inout) :: max_itr1, max_itr2
#endif

    ! work
    real(RP), intent(out) :: TRLP(UKA)    ! Layer temperature at previous step  [K]
    real(RP), intent(out) :: TBLP(UKA)    ! Layer temperature at previous step  [K]
    real(RP), intent(out) :: TGLP(UKA)    ! Layer temperature at previous step  [K]
    real(RP), intent(out) :: A(UKA)
    real(RP), intent(out) :: B(UKA)
    real(RP), intent(out) :: C(UKA)
    real(RP), intent(out) :: D(UKA)
    real(RP), intent(out) :: P(UKA)
    real(RP), intent(out) :: Q(UKA)

    !-- parameters
!    real(RP), parameter     :: SRATIO    = 0.75_RP     ! ratio between direct/total solar [-]
!    real(RP), parameter     :: TFa       = 0.5_RP      ! factor a in Tomita (2009)
!    real(RP), parameter     :: TFb       = 1.1_RP      ! factor b in Tomita (2009)
    real(RP), parameter     :: redf_min  = 1.0E-2_RP   ! minimum reduced factor
    real(RP), parameter     :: redf_max  = 1.0_RP      ! maximum reduced factor
!    real(RP), parameter     :: CAP_water = 4.185E6_RP ! Heat capacity of water (15 deg) [J m-3 K]
!    real(RP), parameter     :: AKS_water = 0.59_RP    ! Thermal conductivity of water   [W m-1 K]

    real(RP), parameter     :: rain_rate_R = 1.0_RP
    real(RP), parameter     :: rain_rate_B = 0.1_RP
    real(RP), parameter     :: rain_rate_G = 0.9_RP

    integer, parameter :: itr_max = 100

    !-- Local variables
!    logical  :: SHADOW = .false.
           !  true  = consider svf and shadow effects,
           !  false = consider svf effect only

    real(RP) :: W, VFGS, VFGW, VFWG, VFWS, VFWW
    real(RP) :: rflux_SW, rflux_LW

    real(RP) :: TRP              ! TRP: at previous time step [K]
    real(RP) :: TBP              ! TBP: at previous time step [K]
    real(RP) :: TGP              ! TGP: at previous time step [K]
    real(RP) :: TCP              ! TCP: at previous time step [K]
    real(RP) :: QCP              ! QCP: at previous time step [kg/kg]
    !
    real(RP) :: RAINRP ! at previous step, rain amount in storage on roof     [kg/m2]
    real(RP) :: RAINBP ! at previous step, rain amount in storage on building [kg/m2]
    real(RP) :: RAINGP ! at previous step, rain amount in storage on road     [kg/m2]

    real(RP) :: RAINT
    real(RP) :: ROFFR, ROFFB, ROFFG ! runoff [kg/m2]

    real(RP) :: LUP, LDN, RUP
    real(RP) :: SUP, SDN

    real(RP) :: LNET, SNET, FLXUV
    real(RP) :: SW                  ! shortwave radition          [W/m/m]
    real(RP) :: LW                  ! longwave radition           [W/m/m]
    real(RP) :: psim,psim2,psim10   ! similality stability shear function for momentum
    real(RP) :: psih,psih2,psih10   ! similality stability shear function for heat

    ! for shadow effect model
    ! real(RP) :: HOUI1, HOUI2, HOUI3, HOUI4, HOUI5, HOUI6, HOUI7, HOUI8
    ! real(RP) :: SLX, SLX1, SLX2, SLX3, SLX4, SLX5, SLX6, SLX7, SLX8
    ! real(RP) :: THEATAZ    ! Solar Zenith Angle [rad]
    ! real(RP) :: THEATAS    ! = PI/2. - THETAZ
    ! real(RP) :: FAI        ! Latitude [rad]

    real(RP) :: SR, SB, SG, RR, RB, RG
    real(RP) :: SR1, SB1, SB2, SG1, SG2
    real(RP) :: RB1, RB2, RG1, RG2
    real(RP) :: HR, EVPR, ELER, G0R, BETR
    real(RP) :: HB, EVPB, ELEB, G0B, BETB
    real(RP) :: HG, EVPG, ELEG, G0G, BETG
    real(RP) :: EVPRp, EVPBp, EVPGp

    real(RP) :: Z
    real(RP) :: QS0R, QS0B, QS0G

    real(RP) :: RIBR, BHR, CDR
    real(RP) :: RIBC, BHC, CDC
    real(RP) :: ALPHAB, ALPHAG, ALPHAC
    real(RP) :: CHR, CHB, CHG, CHC
    real(RP) :: TC1, TC2, QC1, QC2
!    real(RP) :: CAPL1, AKSL1

    real(RP) :: XXX, XXX2, XXX10
    real(RP) :: XXXR ! Monin-Obkhov length for roof [-]
    real(RP) :: XXXC ! Monin-Obkhov length for canopy [-]

    real(RP) :: THA,THC,THS,THS1,THS2
    real(RP) :: RovCP
    real(RP) :: EXN  ! exner function at the surface

    real(RP) :: FracU10, FracT2

    real(RP) :: dt_RP

    integer  :: iteration
    real(RP) :: DTS_MAX_onestep ! DTS_MAX * dt
    real(RP) :: resi1, resi2, resi3 ! residual
    real(RP) :: fact1, fact2, fact3
    real(RP) :: threshold

    ! for Newton method
    real(RP) :: dTR, dTB, dTG, dAC
    real(RP) :: dTRp, dTBp, dTGp, dACp
    real(RP) :: dr1dTR
    real(RP) :: dr1dTB, dr1dTG, dr1dAC, dr2dTB, dr2dTG, dr2dAC, dr3dTB, dr3dTG, dr3dAC
    real(RP) :: dTRdG0R, dTBdG0B, dTGdG0G
    real(RP) :: dG0RdTR, dG0BdTB, dG0BdTG, dG0BdAC, dG0GdTB, dG0GdTG, dG0GdAC
    real(RP) :: dHR, dHBdTB, dHBdTG, dHBdAC, dHGdTB, dHGdTG, dHGdAC
    real(RP) :: dELER, dELEBdTB, dELEBdTG, dELEBdAC, dELEGdTB, dELEGdTG, dELEGdAC
    real(RP) :: dRR
    real(RP) :: dBETR, dBETB, dBETG, BETRP, BETBP, BETGP
    real(RP) :: dRB1dTB, dRB1dTG, dRB2dTB, dRB2dTG
    real(RP) :: dRG1dTB, dRG1dTG, dRG2dTB, dRG2dTG
    real(RP) :: dQCdTB, dQCdTG, dQCdAC
    real(RP) :: dTHCdTB, dTHCdTG, dTHCdAC
    real(RP) :: dALPHACdTB, dALPHACdTG, dALPHACdAC
    real(RP) :: dCHCdTB, dCHCdTG, dCHCdAC
    real(RP) :: CHC_TB, CHC_TG, CHC_AC
    real(RP) :: dQS0R, dQS0B, dQS0G
    real(RP) :: Tdiff, Adiff
    real(RP) :: XXXtmp
    real(RP) :: ALPHACp
    real(RP) :: rdet
    real(RP) :: b1, b2
    real(RP) :: tmp

    integer :: k

    !-----------------------------------------------------------
    ! Set parameters
    !-----------------------------------------------------------

    dt_RP = dt

    threshold = sqrt(EPS)

    RovCP = Rdry / CPdry
    THA   = TA * ( PRE00 / PRSA )**RovCP

    !--- Renew surface and layer temperatures

    TRP = TR
    TBP = TB
    TGP = TG
    TCP = TC
    QCP = QC
    !
    do k = UKS, UKE
       TRLP(k) = TRL(k)
       TBLP(k) = TBL(k)
       TGLP(k) = TGL(k)
    end do
    !


    !--- limiter for surface temp change
    DTS_MAX_onestep = DTS_MAX * dt_RP

    ! "2.0m" has no special meaning, but it is related with BB formulation from Inoue (1963).
    ! Please see subroutine "canopy_wind".
    ! The canopy model is modeled under an assumption that urban canopy lies
    ! below the lowest level of atmospheric model.
    if ( ZDC + Z0C + 2.0_RP >= ZA ) then
       LOG_ERROR("URBAN_DYN_kusaka01_SLC_main",*) 'ZDC + Z0C + 2m must be less than the 1st level! STOP.'
       converged = .false.
#ifdef _OPENACC
       return
#else
       call PRC_abort
#endif
    endif

    W    = 2.0_RP * 1.0_RP * HGT

    rflux_SW   = SSG(1) + SSG(2) ! downward shortwave radiation [W/m2]
    rflux_LW   = LLG(1) + LLG(2) ! downward longwave  radiation [W/m2]

    !--- calculate canopy wind

    call canopy_wind(ZA, UA, Z0C, ZDC, UC)

    !-----------------------------------------------------------
    ! calculate water content (temporally)
    !-----------------------------------------------------------

    RAINT = RAIN * dt_RP ! [kg/m2/s -> kg/m2]

    RAINR = RAINR + RAINT * rain_rate_R
    RAINB = RAINB + RAINT * rain_rate_B
    RAING = RAING + RAINT * rain_rate_G

    VFGS = SVF
    VFGW = 1.0_RP - SVF
    VFWG = ( 1.0_RP - SVF ) * ( 1.0_RP - R ) / W
    VFWS = VFWG
    VFWW = 1.0_RP - 2.0_RP * VFWG

    ! save the initial value
    RAINRP = RAINR
    RAINBP = RAINB
    RAINGP = RAING

    !-----------------------------------------------------------
    ! Radiation : Net Short Wave Radiation at roof/wall/road
    !-----------------------------------------------------------

    if( rflux_SW > 0.0_RP ) then !  SSG is downward short

      ! currently we use no shadow effect model
      !!     IF(.NOT.SHADOW) THEN              ! no shadow effects model

      SR1  = rflux_SW * ( 1.0_RP - ALBR )
      SG1  = rflux_SW * VFGS * ( 1.0_RP - ALBG )
      SB1  = rflux_SW * VFWS * ( 1.0_RP - ALBB )
      SG2  = SB1 * ALBB / ( 1.0_RP - ALBB ) * VFGW * ( 1.0_RP - ALBG )
      SB2  = SG1 * ALBG / ( 1.0_RP - ALBG ) * VFWG * ( 1.0_RP - ALBB )

      SR   = SR1
      SG   = SG1 + SG2
      SB   = SB1 + SB2
      SNET = R * SR + W * SB + RW * SG

    else

      SR   = 0.0_RP
      SG   = 0.0_RP
      SB   = 0.0_RP
      SNET = 0.0_RP

    end if


    EXN = ( PRSS / PRE00 )**RovCP ! exner function

    !-----------------------------------------------------------
    ! Energy balance on roof/wall/road surface
    !-----------------------------------------------------------

    !--------------------------------------------------
    !   Roof
    !--------------------------------------------------

    ! new scheme

    Z   = ZA - ZDC
    BHR = LOG(Z0R/Z0HR) / 0.4_RP

    b1 = CAPR * DZR(1) / dt_RP + 2.0_RP * AKSR / ( DZR(1)+DZR(2) )
    b2 = CAPR * DZR(2) / dt_RP + 2.0_RP * AKSR / ( DZR(1)+DZR(2) ) + 2.0_RP * AKSR / ( DZR(2)+DZR(3) )
    dTRdG0R = 1.0_RP / ( b1 - ( 2.0_RP * AKSR / ( DZR(1)+DZR(2) ) )**2 / b2 )
    ! consider only change at the layers of k<=2
    dTRdG0R = dTRdG0R * 0.5_RP

    EVPR = 0.0
    BETRP = 0.0_RP
    XXXR = 0.0_RP
    dTRp = 1.0e10_RP
    fact1 = 1.0_RP
    !$acc loop seq
    do iteration = 1, itr_max

       THS = TR / EXN ! potential temp

       RIBR = ( GRAV * 2.0_RP / (THA+THS) ) * (THA-THS) * (Z+Z0R) / (UA*UA+EPS)
       call mos(XXXR,CHR,CDR,BHR,RIBR,Z,Z0R,UA,THA,THS,RHOO,i,j)
       ! ignore differential of CHR for Newton method

       call dqs_dtem( TR, RHOS,   & ! [IN]
                      dQS0R, QS0R ) ! [OUT]

       RAINR = max( RAINRP - EVPR * dt_RP, 0.0_RP )
       call cal_beta(BETR, BETR_CONST, RAINR, STRGR)
       dBETR = ( BETR - BETRP ) / sign( max(abs(dTRp), 1E-10_RP), dTRp )
       BETRP = BETR

       RR  = EPSR * ( rflux_LW - STB * (TR**4)  )
       dRR = - EPSR * STB * 4.0_RP * TR**3

       !HR  = RHOO * CPdry * CHR * UA * (TR-TA)
       HR  = RHOO * CPdry * CHR * UA * ( THS - THA ) * EXN
       dHR = RHOO * CPdry * CHR * UA

       EVPR = min( RHOO * CHR * UA * BETR * (QS0R-QA), RAINR / dt_RP )
       ELER = EVPR * LHV
       dELER = RHOO * CHR * UA * ( BETR * dQS0R + dBETR * QS0R ) * LHV

       G0R = SR + RR - HR - ELER
       dG0RdTR = dRR - dHR - dELER

       !--- calculate temperature in roof
       !  if ( STRGR /= 0.0_RP ) then
       !    CAPL1 = CAP_water * (RAINR / (DZR(1) + RAINR)) + CAPR * (DZR(1) / (DZR(1) + RAINR))
       !    AKSL1 = AKS_water * (RAINR / (DZR(1) + RAINR)) + AKSR * (DZR(1) / (DZR(1) + RAINR))
       !  else
       !    CAPL1 = CAPR
       !    AKSL1 = AKSR
       !  endif
       !! 1st layer's cap, aks are replaced.
       !! call multi_layer2(UKE,BOUND,G0R,CAPR,AKSR,TRL,DZR,dt_RP,TRLEND,CAPL1,AKSL1)

       do k = UKS, UKE
          TRL(k) = TRLP(k)
       end do
       call multi_layer(UKE,BOUND, &
            A, B, C, D, P, Q, &
            G0R,CAPR,AKSR,TRL,DZR,dt_RP,TRLEND)
       resi1  = TRL(1) - TR

       if( abs(resi1) < threshold ) then
          TR = TRL(1)
          TR = max( TRP - DTS_MAX_onestep, min( TRP + DTS_MAX_onestep, TR ) )
          exit
       endif

       ! Newton method
       dr1dTR = dTRdG0R * dG0RdTR - 1.0_RP

       dTR = - resi1 / dr1dTR
       dTR = sign( min(abs(dTR), 5.0_RP), dTR )

       if ( iteration > 50 ) then
          if ( dTR*dTRp < 0.0_RP ) then
             fact1 = max( fact1 * 0.5_RP, 0.01_RP )
          else
             fact1 = min( fact1 * 1.5_RP, 1.1_RP )
          end if
       end if
       tmp = TR
       TR = max( TR + dTR * fact1, 100.0_RP )
       dTRp = TR - tmp

       EVPR = ( EVPRp + EVPR ) * 0.5_RP
       EVPRp = EVPR
    enddo

#ifdef DEBUG_KUSAKA01
    cnt_num1 = cnt_num1 + 1
    cnt_itr1 = cnt_itr1 + iteration - 1
    max_itr1 = max(max_itr1, iteration - 1)
#endif

    ! output for debug
    if ( iteration > itr_max .and. debug ) then
       LOG_WARN("URBAN_DYN_kusaka01_SLC_main",*) 'iteration for TR was not converged',PRC_myrank,i,j
       LOG_WARN_CONT(*) '---------------------------------------------------------------------------------'
       LOG_WARN_CONT(*) 'DEBUG Message --- Residual                                          [K] :', resi1
       LOG_WARN_CONT(*) 'DEBUG Message --- TRP : Initial TR                                  [K] :', TRP
#ifdef _OPENACC
       LOG_WARN_CONT(*) 'DEBUG Message --- TRLP: Initial TRL                                 [K] :', TRLP(UKS)
#else
       LOG_WARN_CONT(*) 'DEBUG Message --- TRLP: Initial TRL                                 [K] :', TRLP(:)
#endif
       LOG_WARN_CONT(*) 'DEBUG Message --- rflux_SW  : Shortwave radiation                [W/m2] :', rflux_SW
       LOG_WARN_CONT(*) 'DEBUG Message --- rflux_LW  : Longwave radiation                 [W/m2] :', rflux_LW
       LOG_WARN_CONT(*) 'DEBUG Message --- PRSS: Surface pressure                           [Pa] :', PRSS
       LOG_WARN_CONT(*) 'DEBUG Message --- PRSA: Pressure at 1st atmos layer                 [m] :', PRSA
       LOG_WARN_CONT(*) 'DEBUG Message --- RHOO: Air density                             [kg/m3] :', RHOO
       LOG_WARN_CONT(*) 'DEBUG Message --- RHOS: Surface density                         [kg/m3] :', RHOS
       LOG_WARN_CONT(*) 'DEBUG Message --- RAINRP: Initial RAINR                         [kg/m2] :', RAINRP
       LOG_WARN_CONT(*) 'DEBUG Message --- ZA  : Height at 1st atmos layer                   [m] :', ZA
       LOG_WARN_CONT(*) 'DEBUG Message --- TA  : Temperature at 1st atmos layer              [K] :', TA
       LOG_WARN_CONT(*) 'DEBUG Message --- UA  : Wind speed at 1st atmos layer             [m/s] :', UA
       LOG_WARN_CONT(*) 'DEBUG Message --- QA  : Specific humidity at 1st atmos layer    [kg/kg] :', QA
#ifdef _OPENACC
       LOG_WARN_CONT(*) 'DEBUG Message --- DZR : Depth of surface layer                      [m] :', DZR(1)
#else
       LOG_WARN_CONT(*) 'DEBUG Message --- DZR : Depth of surface layer                      [m] :', DZR(:)
#endif
       LOG_WARN_CONT(*) 'DEBUG Message --- R, W, RW : Normalized height and road width       [-] :', R, W,RW
       LOG_WARN_CONT(*) 'DEBUG Message --- SVF : Sky View Factors                            [-] :', SVF
       LOG_WARN_CONT(*) 'DEBUG Message --- BETR: Evaporation efficiency                      [-] :', BETR
       LOG_WARN_CONT(*) 'DEBUG Message --- EPSR: Surface emissivity of roof                  [-] :', EPSR
       LOG_WARN_CONT(*) 'DEBUG Message --- CAPR: Heat capacity of roof                 [J m-3 K] :', CAPR
       LOG_WARN_CONT(*) 'DEBUG Message --- AKSR: Thermal conductivity of roof          [W m-1 K] :', AKSR
       LOG_WARN_CONT(*) 'DEBUG Message --- QS0R: Surface specific humidity               [kg/kg] :', QS0R
       LOG_WARN_CONT(*) 'DEBUG Message --- ZDC : Desplacement height of canopy               [m] :', ZDC
       LOG_WARN_CONT(*) 'DEBUG Message --- Z0R : Momentum roughness length of roof           [m] :', Z0R
       LOG_WARN_CONT(*) 'DEBUG Message --- Z0HR: Thermal roughness length of roof            [m] :', Z0HR
       LOG_WARN_CONT(*) '---------------------------------------------------------------------------------'
     endif

    !--- update only fluxes ----
     THS   = TR / EXN
     RIBR = ( GRAV * 2.0_RP / (THA+THS) ) * (THA-THS) * (Z+Z0R) / (UA*UA+EPS)
     call mos(XXXR,CHR,CDR,BHR,RIBR,Z,Z0R,UA,THA,THS,RHOO,i,j)

     call qsat( TR, RHOS, & ! [IN]
                QS0R      ) ! [OUT]

     call cal_beta(BETR, BETR_CONST, RAINR, STRGR)

     RR      = EPSR * ( rflux_LW - STB * (TR**4) )
     HR      = RHOO * CPdry * CHR * UA * (THS-THA) * EXN
     EVPR    = min( RHOO * CHR * UA * BETR * (QS0R-QA), RAINR / dt_RP )
     ELER    = EVPR * LHV

!     G0R     = SR + RR - HR - ELER + EFLX
     G0R     = SR + RR - HR - ELER
     RAINR   = max( RAINRP - EVPR * dt_RP, 0.0_RP )

     do k = UKS, UKE
        TRL(k)  = TRLP(k)
     end do
     call multi_layer(UKE,BOUND, &
          A, B, C, D, P, Q, &
          G0R,CAPR,AKSR,TRL,DZR,dt_RP,TRLEND)
     resi1   = TRL(1) - TR
     TR      = TRL(1)

     if ( abs(resi1) > DTS_MAX_onestep ) then
       if ( abs(resi1) > DTS_MAX_onestep*10.0_RP ) then
         LOG_ERROR("URBAN_DYN_Kusaka01_main",*) 'tendency of TR exceeded a limit! STOP.'
         LOG_ERROR_CONT(*) 'previous TR and updated TR(TRL(1)) is ',TR-resi1, TR
         converged = .false.
#ifdef _OPENACC
         return
#else
         call PRC_abort
#endif
       endif
       LOG_WARN("URBAN_DYN_Kusaka01_main",*) 'tendency of TR exceeded a limit'
       LOG_WARN_CONT(*) 'previous TR and updated TR(TRL(1)) is ', TR-resi1, TR
     endif

    !--------------------------------------------------
    !   Wall and Road
    !--------------------------------------------------

    ! new scheme

    ! empirical form
     ALPHAB = 6.15_RP + 4.18_RP * UC
     if( UC > 5.0_RP ) ALPHAB = 7.51_RP * (UC**0.78_RP )
     ALPHAG = 6.15_RP + 4.18_RP * UC
     if( UC > 5.0_RP ) ALPHAG = 7.51_RP * (UC**0.78_RP )
     CHB = ALPHAB / RHOO / CPdry / UC
     CHG = ALPHAG / RHOO / CPdry / UC

     Z   = ZA - ZDC
     THC = TC / EXN
     BHC = LOG(Z0C/Z0HC) / 0.4_RP

     b1 = CAPB * DZB(1) / dt_RP + 2.0_RP * AKSB / ( DZB(1)+DZB(2) )
     b2 = CAPB * DZB(2) / dt_RP + 2.0_RP * AKSB / ( DZB(1)+DZB(2) ) + 2.0_RP * AKSB / ( DZB(2)+DZB(3) )
     dTBdG0B = 1.0_RP / ( b1 - ( 2.0_RP * AKSB / ( DZB(1)+DZB(2) ) )**2 / b2 )
     ! consider only change at the layers of k<=2
     dTBdG0B = dTBdG0B * 0.5_RP

     b1 = CAPG * DZG(1) / dt_RP + 2.0_RP * AKSG / ( DZG(1)+DZG(2) )
     b2 = CAPG * DZG(2) / dt_RP + 2.0_RP * AKSG / ( DZG(1)+DZG(2) ) + 2.0_RP * AKSG / ( DZG(2)+DZG(3) )
     dTGdG0G = 1.0_RP / ( b1 - ( 2.0_RP * AKSG / ( DZG(1)+DZG(2) ) )**2 / b2 )
     ! consider only change at the layers of k<=2
     dTGdG0G = dTGdG0G * 0.5_RP

     EVPB = 0.0_RP
     EVPG = 0.0_RP
     BETBP = 0.0_RP
     BETGP = 0.0_RP
     XXXC = 0.0_RP
     ALPHACp = ( ALPHAB + ALPHAG ) * 0.5_RP
     ALPHAC = ALPHACp
     dTBp = 1.0e10_RP
     dTGp = 1.0e10_RP
     dACp = 1.0_RP
     fact1 = 1.0_RP
     fact2 = 1.0_RP
     fact3 = 1.0_RP
     !$acc loop seq
     do iteration = 1, itr_max

        THS1   = TB / EXN
        THS2   = TG / EXN

        if ( iteration > 1 ) then

           Tdiff = TB * sqrt(EPS) * 2.0_RP
           Adiff = sign( ALPHAC * sqrt(EPS) * 2.0_RP, dAC )

           ! TB = TB + Tdiff
           TC1    =  RW*ALPHAC    + RW*ALPHAG    + W*ALPHAB
           TC2    =  RW*ALPHAC*THA + W*ALPHAB*(TB+Tdiff)/EXN + RW*ALPHAG*THS2
           THC    =  TC2 / TC1
           RIBC = ( GRAV * 2.0_RP / (THA+THC) ) * (THA-THC) * (Z+Z0C) / (UA*UA+EPS)
           XXXtmp = XXXC
           RIBC = ( GRAV * 2.0_RP / (THA+THC) ) * (THA-THC) * (Z+Z0C) / (UA*UA+EPS)
           call mos(XXXtmp,CHC_TB,CDC,BHC,RIBC,Z,Z0C,UA,THA,THC,RHOO,i,j)

           ! TG = TG + Tdiff
           TC1    =  RW*ALPHAC    + RW*ALPHAG    + W*ALPHAB
           TC2    =  RW*ALPHAC*THA + W*ALPHAB*THS1 + RW*ALPHAG*(TG+Tdiff)/EXN
           THC    =  TC2 / TC1
           RIBC = ( GRAV * 2.0_RP / (THA+THC) ) * (THA-THC) * (Z+Z0C) / (UA*UA+EPS)
           XXXtmp = XXXC
           call mos(XXXtmp,CHC_TG,CDC,BHC,RIBC,Z,Z0C,UA,THA,THC,RHOO,i,j)

           ! ALPHAC = ALPHAC + Adiff
           TC1    =  RW*(ALPHAC+Adiff)    + RW*ALPHAG    + W*ALPHAB
           TC2    =  RW*(ALPHAC+Adiff)*THA + W*ALPHAB*THS1 + RW*ALPHAG*THS2
           THC    =  TC2 / TC1
           RIBC = ( GRAV * 2.0_RP / (THA+THC) ) * (THA-THC) * (Z+Z0C) / (UA*UA+EPS)
           XXXtmp = XXXC
           call mos(XXXtmp,CHC_AC,CDC,BHC,RIBC,Z,Z0C,UA,THA,THC,RHOO,i,j)

        end if

        TC1    =  RW*ALPHAC    + RW*ALPHAG    + W*ALPHAB
        !TC2    =  RW*ALPHAC*THA + RW*ALPHAG*TG + W*ALPHAB*TB
        TC2    =  RW*ALPHAC*THA + W*ALPHAB*THS1 + RW*ALPHAG*THS2
        THC    =  TC2 / TC1
        RIBC = ( GRAV * 2.0_RP / (THA+THC) ) * (THA-THC) * (Z+Z0C) / (UA*UA+EPS)
        call mos(XXXC,CHC,CDC,BHC,RIBC,Z,Z0C,UA,THA,THC,RHOO,i,j)

        if ( iteration > 1 ) then
           dCHCdTB = ( CHC_TB - CHC ) / Tdiff
           dCHCdTG = ( CHC_TG - CHC ) / Tdiff
           dCHCdAC = ( CHC_AC - CHC ) / Adiff
        else
           dCHCdTB = 0.0_RP
           dCHCdTG = 0.0_RP
           dCHCdAC = 0.0_RP
        end if

        ALPHAC = CHC * RHOO * CPdry * UA
        dALPHACdTB = dCHCdTB * RHOO * CPdry * UA
        dALPHACdTG = dCHCdTG * RHOO * CPdry * UA
        dALPHACdAC = dCHCdAC * RHOO * CPdry * UA

        call dqs_dtem( TB, RHOS, dQS0B, QS0B )
        call dqs_dtem( TG, RHOS, dQS0G, QS0G )

        RAINB = max( ( RAINBP - EVPB * dt_RP ), 0.0_RP )
        RAING = max( ( RAINGP - EVPG * dt_RP ), 0.0_RP )
        call cal_beta(BETB, BETB_CONST, RAINB, STRGB)
        call cal_beta(BETG, BETG_CONST, RAING, STRGG)
        dBETB = ( BETB - BETBP ) / sign( max(abs(dTBp),1E-10_RP), dTBp )
        dBETG = ( BETG - BETGP ) / sign( max(abs(dTGp),1E-10_RP), dTGp )
        BETBP = BETB
        BETGP = BETG


        TC1 = RW*ALPHAC    + RW*ALPHAG    + W*ALPHAB
        !TC2 = RW*ALPHAC*TA + RW*ALPHAG*TG + W*ALPHAB*TB
        TC2 = RW*ALPHAC*THA + W*ALPHAB*THS1 + RW*ALPHAG*THS2
        THC = TC2 / TC1
        dTHCdTB = ( ( RW*dALPHACdTB*THA +  W*ALPHAB/EXN ) * TC1 - RW*dALPHACdTB * TC2 ) / TC1**2
        dTHCdTG = ( ( RW*dALPHACdTG*THA + RW*ALPHAG/EXN ) * TC1 - RW*dALPHACdTG * TC2 ) / TC1**2
        dTHCdAC = ( ( RW*           THA                 ) * TC1 - RW            * TC2 ) / TC1**2

        QC1 = RW*(CHC*UA)    + RW*(CHG*BETG*UC)      + W*(CHB*BETB*UC)
        QC2 = RW*(CHC*UA)*QA + RW*(CHG*BETG*UC)*QS0G + W*(CHB*BETB*UC)*QS0B
        QC  = max( QC2 / ( QC1 + EPS ), 0.0_RP )
        dQCdTB = ( ( RW*(dCHCdTB*UA)*QA +  W*(CHB*UC)*(BETB*dQS0B+dBETB*QS0B) ) * QC1 &
                 - ( RW*(dCHCdTB*UA)    +  W*(CHB*UC*dBETB) ) * QC2 ) / ( QC1**2 + EPS )
        dQCdTG = ( ( RW*(dCHCdTG*UA)*QA + RW*(CHG*UC)*(BETG*dQS0G+dBETG*QS0G) ) * QC1 &
                 - ( RW*(dCHCdTG*UA)    + RW*(CHG*UC*dBETG) ) * QC2 ) / ( QC1**2 + EPS )
        dQCdAC = ( ( RW*(dCHCdAC*UA)*QA                                       ) * QC1 &
                 - ( RW*(dCHCdAC*UA)                        ) * QC2 ) / ( QC1**2 + EPS )

        RG1 = EPSG * ( rflux_LW * VFGS            &
                     + EPSB * VFGW * STB * TB**4  &
                     - STB * TG**4                )
        dRG1dTB = EPSG * EPSB * VFGW * STB * 4.0_RP * TB**3
        dRG1dTG = - EPSG * STB * 4.0_RP * TG**3

        RB1 = EPSB * ( rflux_LW * VFWS            &
                     + EPSG * VFWG * STB * TG**4  &
                     + EPSB * VFWW * STB * TB**4  &
                     - STB * TB**4                )
        dRB1dTB = EPSB * ( EPSB * VFWW - 1.0_RP ) * STB * 4.0_RP * TB**3
        dRB1dTG = EPSB * EPSG * VFWG * STB * 4.0_RP * TG**3

        RG2 = EPSG * ( (1.0_RP-EPSB) * VFGW * VFWS * rflux_LW            &
                     + (1.0_RP-EPSB) * VFGW * VFWG * EPSG * STB * TG**4  &
                     + EPSB * (1.0_RP-EPSB) * VFGW * VFWW * STB * TB**4  )
        dRG2dTB = EPSG * EPSB * (1.0_RP-EPSB) * VFGW * VFWW * STB * 4.0_RP * TB**3
        dRG2dTG = EPSG * (1.0_RP-EPSB) * VFGW * VFWG * EPSG * STB * 4.0_RP * TG**3

        RB2 = EPSB * ( (1.0_RP-EPSG) * VFWG * VFGS * rflux_LW                            &
                     + (1.0_RP-EPSB) * VFWS * VFWW * rflux_LW                            &
                     + (1.0_RP-EPSB) * EPSG * VFWG * VFWW * STB * TG**4                  &
                     + (1.0_RP-EPSG) * EPSB * VFGW * VFWG * STB * TB**4                  &
                     + (1.0_RP-EPSB) * EPSB * VFWW * (1.0_RP-2.0_RP*VFWS) * STB * TB**4  )
        dRB2dTB = EPSB**2 * ( (1.0_RP-EPSG) * VFGW * VFWG + (1.0_RP-EPSB) * VFWW * (1.0_RP-2.0_RP*VFWS) ) * STB * 4.0_RP * TB**3
        dRB2dTG = EPSB * (1.0_RP-EPSB) * EPSG * VFWG * VFWW * STB * 4.0_RP * TG**3

        RG = RG1 + RG2
        RB = RB1 + RB2

        HB     = RHOO * CPdry * CHB * UC * ( THS1 - THC ) * EXN
        dHBdTB = RHOO * CPdry * CHB * UC * ( 1.0_RP - dTHCdTB * EXN )
        dHBdTG = RHOO * CPdry * CHB * UC * (        - dTHCdTG * EXN )
        dHBdAC = RHOO * CPdry * CHB * UC * (        - dTHCdAC * EXN )

        EVPB = min( RHOO * CHB * UC * BETB * (QS0B-QC), RAINB / dt_RP )
        ELEB = EVPB * LHV
        dELEBdTB = RHOO * CHB * UC * ( BETB * ( dQS0B - dQCdTB ) + dBETB * (QS0B-QC) ) * LHV
        dELEBdTG = RHOO * CHB * UC *   BETB * (       - dQCdTG ) * LHV
        dELEBdAC = RHOO * CHB * UC *   BETB * (       - dQCdAC ) * LHV

!        G0B = SB + RB - HB - ELEB + EFLX
        G0B = SB + RB - HB - ELEB
        dG0BdTB = dRB1dTB + dRB2dTB - dHBdTB - dELEBdTB
        dG0BdTG = dRB1dTG + dRB2dTG - dHBdTG - dELEBdTG
        dG0BdAC =                   - dHBdAC - dELEBdAC

        HG     = RHOO * CPdry * CHG * UC * ( THS2 - THC ) * EXN
        dHGdTB = RHOO * CPdry * CHG * UC * (        - dTHCdTB * EXN )
        dHGdTG = RHOO * CPdry * CHG * UC * ( 1.0_RP - dTHCdTG * EXN )
        dHGdAC = RHOO * CPdry * CHG * UC * (        - dTHCdAC * EXN )

        EVPG = min( RHOO * CHG * UC * BETG * (QS0G-QC), RAING / dt_RP )
        ELEG = EVPG * LHV
        dELEGdTB = RHOO * CHG * UC *   BETG * (       - dQCdTB ) * LHV
        dELEGdTG = RHOO * CHG * UC * ( BETG * ( dQS0G - dQCdTG ) + dBETG * (QS0G-QC) ) * LHV
        dELEGdAC = RHOO * CHG * UC *   BETG * (       - dQCdAC ) * LHV

!        G0G = SG + RG - HG - ELEG + EFLX
        G0G = SG + RG - HG - ELEG
        dG0GdTB = dRG1dTB + dRG2dTB - dHGdTB - dELEGdTB
        dG0GdTG = dRG1dTG + dRG2dTG - dHGdTG - dELEGdTG
        dG0GdAC =                   - dHGdAC - dELEGdAC

        do k = UKS, UKE
           TBL(k) = TBLP(k)
        end do
        call multi_layer(UKE,BOUND, &
           A, B, C, D, P, Q, &
           G0B,CAPB,AKSB,TBL,DZB,dt_RP,TBLEND)
        resi1  = TBL(1) - TB

        do k = UKS, UKE
           TGL(k) = TGLP(k)
        end do
        call multi_layer(UKE,BOUND, &
           A, B, C, D, P, Q, &
           G0G,CAPG,AKSG,TGL,DZG,dt_RP,TGLEND)
        resi2 = TGL(1) - TG


        resi3 = ALPHAC - ALPHACp

        if ( abs(resi1) < threshold .AND. abs(resi2) < threshold .AND. abs(resi3) < threshold ) then
           TB = TBL(1)
           TG = TGL(1)
           TB = max( TBP - DTS_MAX_onestep, min( TBP + DTS_MAX_onestep, TB ) )
           TG = max( TGP - DTS_MAX_onestep, min( TGP + DTS_MAX_onestep, TG ) )
           exit
        endif

        ! Newton method
        dr1dTB = dTBdG0B * dG0BdTB - 1.0_RP
        dr1dTG = dTBdG0B * dG0BdTG
        dr1dAC = dTBdG0B * dG0BdAC
        dr2dTB = dTGdG0G * dG0GdTB
        dr2dTG = dTGdG0G * dG0GdTG - 1.0_RP
        dr2dAC = dTGdG0G * dG0GdAC
        dr3dTB = dALPHACdTB
        dr3dTG = dALPHACdTG
        dr3dAC = dALPHACdAC - 1.0_RP
        if ( iteration > 3 ) dr3dAC = min( dr3dAC, -0.02_RP )

        rdet = 1.0_RP &
             / ( dr1dTB * dr2dTG * dr3dAC + dr1dTG * dr2dAC * dr3dTB + dr1dAC * dr2dTB * dr3dTG &
               - dr1dAC * dr2dTG * dr3dTB - dr1dTG * dr2dTB * dr3dAC - dr1dTB * dr2dAC * dr3dTG )
        dTB = ( - ( dr2dTG * dr3dAC - dr2dAC * dr3dTG ) * resi1 &
                + ( dr1dTG * dr3dAC - dr1dAC * dr3dTG ) * resi2 &
                - ( dr1dTG * dr2dAC - dr1dAC * dr2dTG ) * resi3 ) * rdet
        dTG = (   ( dr2dTB * dr3dAC - dr2dAC * dr3dTB ) * resi1 &
                - ( dr1dTB * dr3dAC - dr1dAC * dr3dTB ) * resi2 &
                + ( dr1dTB * dr2dAC - dr1dAC * dr2dTB ) * resi3 ) * rdet
        dAC = ( - ( dr2dTB * dr3dTG - dr2dTG * dr3dTB ) * resi1 &
                + ( dr1dTB * dr3dTG - dr1dTG * dr3dTB ) * resi2 &
                - ( dr1dTB * dr2dTG - dr1dTG * dr2dTB ) * resi3 ) * rdet

        dTB = sign( min(abs(dTB), 5.0_RP), dTB )
        dTG = sign( min(abs(dTG), 5.0_RP), dTG )
        dAC = sign( min(abs(dAC), 50.0_RP), dAC )

        if ( iteration > 50 ) then
           if ( dTB*dTBp < 0.0_RP ) then
              fact1 = max( fact1 * 0.5_RP, 1E-10_RP )
           else
              fact1 = min( fact1 * 1.5_RP, 1.1_RP )
           end if
           if ( dTG*dTGp < 0.0_RP ) then
              fact2 = max( fact2 * 0.5_RP, 1E-10_RP )
           else
              fact2 = min( fact2 * 1.5_RP, 1.1_RP )
           end if
           if ( dAC*dACp < 0.0_RP ) then
              fact3 = max( fact3 * 0.5_RP, 1E-10_RP )
           else
              fact3 = min( fact3 * 1.5_RP, 1.1_RP )
           end if
        end if

        tmp = TB
        TB = max( TB + dTB * fact1, 100.0_RP )
        dTBp = TB - tmp
        tmp = TG
        TG = max( TG + dTG * fact2, 100.0_RP )
        dTGp = TG - tmp
        ALPHAC = max( ALPHACp + dAC * fact3, EPS )
        dACp = ALPHAC - ALPHACp

        if ( abs(dTBp) < threshold .AND. abs(dTGp) < threshold .AND. abs(dACp) < threshold ) then
           TB = TBL(1)
           TG = TGL(1)
           TB = max( TBP - DTS_MAX_onestep, min( TBP + DTS_MAX_onestep, TB ) )
           TG = max( TGP - DTS_MAX_onestep, min( TGP + DTS_MAX_onestep, TG ) )
           exit
        endif

        ALPHACp = ALPHAC

        EVPB = ( EVPBp + EVPB ) * 0.5_RP
        EVPBp = EVPB

        EVPG = ( EVPGp + EVPG ) * 0.5_RP
        EVPGp = EVPG

     enddo

#ifdef DEBUG_KUSAKA01
     cnt_num2 = cnt_num2 + 1
     cnt_itr2 = cnt_itr2 + iteration - 1
     max_itr2 = max(max_itr2, iteration - 1)
#endif

     ! output for debug
     if ( iteration > itr_max .and. debug ) then
       LOG_WARN("URBAN_DYN_Kusaka01_main",*) 'iteration for TB/TG was not converged',PRC_myrank,i,j
       LOG_WARN_CONT(*) '---------------------------------------------------------------------------------'
       LOG_WARN_CONT(*) 'DEBUG Message --- Residual                                        [K] :', resi1, resi2, resi3
       LOG_WARN_CONT(*) 'DEBUG Message --- TBP : Initial TB                                [K] :', TBP
#ifdef _OPENACC
       LOG_WARN_CONT(*) 'DEBUG Message --- TBLP: Initial TBL                               [K] :', TBLP(UKS)
#else
       LOG_WARN_CONT(*) 'DEBUG Message --- TBLP: Initial TBL                               [K] :', TBLP(:)
#endif
       LOG_WARN_CONT(*) 'DEBUG Message --- TGP : Initial TG                                [K] :', TGP
#ifdef _OPENACC
       LOG_WARN_CONT(*) 'DEBUG Message --- TGLP: Initial TGL                               [K] :', TGLP(UKS)
#else
       LOG_WARN_CONT(*) 'DEBUG Message --- TGLP: Initial TGL                               [K] :', TGLP(:)
#endif
       LOG_WARN_CONT(*) 'DEBUG Message --- TCP : Initial TC                                [K] :', TCP
       LOG_WARN_CONT(*) 'DEBUG Message --- QCP : Initial QC                                [K] :', QCP
       LOG_WARN_CONT(*) 'DEBUG Message --- UC  : Canopy wind                             [m/s] :', UC
       LOG_WARN_CONT(*) 'DEBUG Message --- rflux_SW  : Shortwave radiation              [W/m2] :', rflux_SW
       LOG_WARN_CONT(*) 'DEBUG Message --- rflux_LW  : Longwave radiation               [W/m2] :', rflux_LW
       LOG_WARN_CONT(*) 'DEBUG Message --- PRSS: Surface pressure                         [Pa] :', PRSS
       LOG_WARN_CONT(*) 'DEBUG Message --- PRSA: Pressure at 1st atmos layer               [m] :', PRSA
       LOG_WARN_CONT(*) 'DEBUG Message --- RHOO: Air density                           [kg/m3] :', RHOO
       LOG_WARN_CONT(*) 'DEBUG Message --- RHOS: Surface density                       [kg/m3] :', RHOS
       LOG_WARN_CONT(*) 'DEBUG Message --- RAINBP: Initial RAINB                       [kg/m2] :', RAINBP
       LOG_WARN_CONT(*) 'DEBUG Message --- RAINGP: Initial RAING                       [kg/m2] :', RAINGP
       LOG_WARN_CONT(*) 'DEBUG Message --- ZA  : Height at 1st atmos layer                 [m] :', ZA
       LOG_WARN_CONT(*) 'DEBUG Message --- TA  : Temperature at 1st atmos layer            [K] :', TA
       LOG_WARN_CONT(*) 'DEBUG Message --- UA  : Wind speed at 1st atmos layer           [m/s] :', UA
       LOG_WARN_CONT(*) 'DEBUG Message --- QA  : Specific humidity at 1st atmos layer  [kg/kg] :', QA
#ifdef _OPENACC
       LOG_WARN_CONT(*) 'DEBUG Message --- DZB : Depth of surface layer                    [m] :', DZB(1)
       LOG_WARN_CONT(*) 'DEBUG Message --- DZG : Depth of surface layer                    [m] :', DZG(1)
#else
       LOG_WARN_CONT(*) 'DEBUG Message --- DZB : Depth of surface layer                    [m] :', DZB(:)
       LOG_WARN_CONT(*) 'DEBUG Message --- DZG : Depth of surface layer                    [m] :', DZG(:)
#endif
       LOG_WARN_CONT(*) 'DEBUG Message --- R, W, RW  : Normalized height and road width    [-] :', R, W,RW
       LOG_WARN_CONT(*) 'DEBUG Message --- SVF       : Sky View Factors                    [-] :', SVF
       LOG_WARN_CONT(*) 'DEBUG Message --- BETB,BETG : Evaporation efficiency              [-] :', BETB,BETG
       LOG_WARN_CONT(*) 'DEBUG Message --- EPSB,EPSG : Surface emissivity                  [-] :', EPSB,EPSG
       LOG_WARN_CONT(*) 'DEBUG Message --- CAPB,CAPG : Heat capacity                 [J m-3 K] :', CAPB,CAPG
       LOG_WARN_CONT(*) 'DEBUG Message --- AKSB,AKSG : Thermal conductivity          [W m-1 K] :', AKSB,AKSB
       LOG_WARN_CONT(*) 'DEBUG Message --- QS0B,QS0G : Surface specific humidity       [kg/kg] :', QS0B,QS0G
       LOG_WARN_CONT(*) 'DEBUG Message --- ZDC       : Desplacement height of canopy       [m] :', ZDC
       LOG_WARN_CONT(*) 'DEBUG Message --- Z0M       : Momentum roughness length of canopy [m] :', Z0C
       LOG_WARN_CONT(*) 'DEBUG Message --- Z0H/Z0E   : Thermal roughness length of canopy  [m] :', Z0HC
       LOG_WARN_CONT(*) 'DEBUG Message --- LHV       : Latent heat of vapor           [J kg-1] :', LHV
       LOG_WARN_CONT(*) 'DEBUG Message --- dt        : Time step                           [s] :', dt_RP
       LOG_WARN_CONT(*) '---------------------------------------------------------------------------------'
     endif

     !--- update only fluxes ----

     ! this is for TC, QC
     call qsat( TB, RHOS, QS0B )
     call qsat( TG, RHOS, QS0G )

     call cal_beta(BETB, BETB_CONST, RAINB, STRGB)
     call cal_beta(BETG, BETG_CONST, RAING, STRGG)


     RG1      = EPSG * ( rflux_LW * VFGS                  &
                       + EPSB * VFGW * STB * TB**4  &
                       - STB * TG**4                )
     RB1      = EPSB * ( rflux_LW * VFWS                  &
                       + EPSG * VFWG * STB * TG**4  &
                       + EPSB * VFWW * STB * TB**4  &
                       - STB * TB**4                )

     RG2      = EPSG * ( (1.0_RP-EPSB) * VFGW * VFWS * rflux_LW                  &
                       + (1.0_RP-EPSB) * VFGW * VFWG * EPSG * STB * TG**4  &
                       + EPSB * (1.0_RP-EPSB) * VFGW * VFWW * STB * TB**4  )
     RB2      = EPSB * ( (1.0_RP-EPSG) * VFWG * VFGS * rflux_LW                                 &
                       + (1.0_RP-EPSG) * EPSB * VFGW * VFWG * STB * TB**4                 &
                       + (1.0_RP-EPSB) * VFWS * VFWW * rflux_LW                                 &
                       + (1.0_RP-EPSB) * VFWG * VFWW * STB * EPSG * TG**4                 &
                       + EPSB * (1.0_RP-EPSB) * VFWW * (1.0_RP-2.0_RP*VFWS) * STB * TB**4 )

     RG       = RG1 + RG2
     RB       = RB1 + RB2

     THS1   = TB / EXN
     THS2   = TG / EXN
     THC    = TC / EXN

     HB    = RHOO * CPdry * CHB * UC * (THS1-THC) * EXN
     EVPB  = min( RHOO * CHB * UC * BETB * (QS0B-QC), RAINB / dt_RP )
     ELEB  = EVPB * LHV
!     G0B   = SB + RB - HB - ELEB + EFLX
     G0B   = SB + RB - HB - ELEB
     RAINB = max( RAINBP - EVPB * dt_RP, 0.0_RP )

     HG    = RHOO * CPdry * CHG * UC * (THS2-THC) * EXN
     EVPG  = min( RHOO * CHG * UC * BETG * (QS0G-QC), RAING / dt_RP )
     ELEG  = EVPG * LHV
!     G0G   = SG + RG - HG - ELEG + EFLX
     G0G   = SG + RG - HG - ELEG
     RAING = max( RAINGP - EVPG * dt_RP, 0.0_RP )

     do k = UKS, UKE
        TBL(k) = TBLP(k)
     end do
     call multi_layer(UKE,BOUND, &
          A, B, C, D, P, Q, &
          G0B,CAPB,AKSB,TBL,DZB,dt_RP,TBLEND)
     resi1  = TBL(1) - TB
     TB     = TBL(1)

     do k = UKS, UKE
        TGL(k) = TGLP(k)
     end do
     call multi_layer(UKE,BOUND, &
          A, B, C, D, P, Q, &
          G0G,CAPG,AKSG,TGL,DZG,dt_RP,TGLEND)
     resi2  = TGL(1) - TG
     TG     = TGL(1)

     if ( abs(resi1) > DTS_MAX_onestep ) then
        if ( abs(resi1) > DTS_MAX_onestep*10.0_RP ) then
           LOG_ERROR("URBAN_DYN_Kusaka01_main",*) 'tendency of TB exceeded a limit! STOP.'
           LOG_ERROR_CONT(*) 'previous TB and updated TB(TBL(1)) is ', TB-resi1,TB
           converged = .false.
#ifdef _OPENACC
           return
#else
           call PRC_abort
#endif
        endif
        LOG_WARN("URBAN_DYN_Kusaka01_main",*) 'tendency of TB exceeded a limit'
        LOG_WARN_CONT(*) 'previous TB and updated TB(TBL(1)) is ', TB-resi1, TB
     endif

     if ( abs(resi2) > DTS_MAX_onestep ) then
        if ( abs(resi2) > DTS_MAX_onestep*10.0_RP ) then
           LOG_ERROR("URBAN_DYN_Kusaka01_main",*) 'tendency of TG exceeded a limit! STOP.'
           LOG_ERROR_CONT(*) 'previous TG and updated TG(TGL(1)) is ', TG-resi2, TG, resi2
           converged = .false.
#ifdef _OPENACC
           return
#else
           call PRC_abort
#endif
        endif
        LOG_WARN("URBAN_DYN_Kusaka01_main",*) 'tendency of TG exceeded a limit'
        LOG_WARN_CONT(*) 'previous TG and updated TG(TGL(1)) is ', TG-resi2, TG
     endif

    !-----------------------------------------------------------
    ! Total Fluxes from Urban Canopy
    !-----------------------------------------------------------

    FLXUV = ( R*CDR + RW*CDC ) * UA * UA
    MFLX  = - RHOO * FLXUV                ! Momentum flux      [kg/m/s2]
    SH    = ( R*HR   + W*HB   + RW*HG )   ! Sensible heat flux [W/m/m]
    LH    = ( R*ELER + W*ELEB + RW*ELEG ) ! Latent heat flux   [W/m/m]
    GHFLX = R*G0R + W*G0B + RW*G0G
    LNET  = R*RR + W*RB + RW*RG

    !-----------------------------------------------------------
    ! Grid average
    !-----------------------------------------------------------

    LW = rflux_LW - LNET     ! Upward longwave radiation   [W/m/m]
    SW = rflux_SW - SNET     ! Upward shortwave radiation  [W/m/m]
    RN = (SNET+LNET)    ! Net radiation [W/m/m]

    !--- shortwave radiation
    SDN =  R + W * (VFWS + VFGS * ALBG * VFWG)  + RW * (VFGS + VFWS * ALBB * VFGW)
    SUP =  R * ALBR                                        &
         + W *  ( VFWS * ALBB + VFGS * ALBG * VFWG *ALBB ) &
         + RW * ( VFGS * ALBG + VFWS * ALBB * VFGW *ALBG )

    ALBD_SW_grid = SUP / SDN

    !--- longwave radiation
    LDN =  R + W*VFWS + RW*VFGS
    LUP =  R * (1.0_RP-EPSR)                                                           &
         + W*( (1.0_RP-EPSB*VFWW)*(1.0_RP-EPSB)*VFWS - EPSB*VFWG*(1.0_RP-EPSG)*VFGS )  &
         + RW*( (1.0_RP-EPSG)*VFGS - EPSG*(1.0_RP-VFGS)*(1.0_RP-EPSB)*VFWS )

    RUP = (LDN - LUP) * rflux_LW - LNET
    ALBD_LW_grid = LUP / LDN


    ! RUP =  R * (EPSR * STB * TR**4 ) &
    !      + W * (EPSB*STB*(TB**4) - EPSB*EPSG*VFWG*STB*(TG**4) - EPSB*EPSB*VFWW*STB*(TB**4)  &
    !           - EPSB *(1.0_RP-EPSG) * EPSB * VFGW * VFWG * STB * (TB**4)         &
    !           - EPSB *(1.0_RP-EPSB) * VFWG * VFWW * STB * EPSG * (TG**4)         &
    !           - EPSB * EPSB * (1.0_RP-EPSB) * VFWW * VFWW * STB * (TB**4) )      &
    !      + RW * (EPSG*STB*(TG**4) - EPSG * EPSB * VFGW * STB * (TB**4)             &
    !           - EPSG * EPSB * (1.0_RP-EPSB) * (1.0_RP-SVF) * VFWW * STB * TB**4  )


    SHR = HR             ! Sensible heat flux on roof [W/m/m]
    SHB = HB             ! Sensible heat flux on wall [W/m/m]
    SHG = HG             ! Sensible heat flux on road [W/m/m]
    LHR = ELER           ! Latent heat flux on road [W/m/m]
    LHB = ELEB           ! Latent heat flux on wall [W/m/m]
    LHG = ELEG           ! Latent heat flux on road [W/m/m]
    GHR = G0R            ! Ground heat flux on roof [W/m/m]
    GHB = G0B            ! Ground heat flux on wall [W/m/m]
    GHG = G0G            ! Ground heat flux on road [W/m/m]
    RNR = SR + RR        ! Net radiation on roof [W/m/m]
    RNB = SB + RB        ! Net radiation on building [W/m/m]
    RNG = SG + RG        ! Net radiation on ground [W/m/m]

    !-----------------------------------------------------------
    !  calculate the ranoff
    !-----------------------------------------------------------
    ROFFR = max( RAINR - STRGR, 0.0_RP )
    ROFFB = max( RAINB - STRGB, 0.0_RP )
    ROFFG = max( RAING - STRGG, 0.0_RP )

    ROFF = R * ROFFR + RW * ( ROFFB + ROFFG )

    RAINR = max( RAINR - ROFFR, 0.0_RP )
    RAINB = max( RAINB - ROFFB, 0.0_RP )
    RAING = max( RAING - ROFFG, 0.0_RP )

    !-----------------------------------------------------------
    !  diagnostic GRID AVERAGED TS from upward logwave
    !-----------------------------------------------------------

    RTS = ( RUP / STB / ( 1.0_RP-ALBD_LW_grid) )**0.25

    !-----------------------------------------------------------
    !  diagnostic grid average U10, V10, T2, Q2 from urban
    !  Below method would be better to be improved. This is tentative method.
    !-----------------------------------------------------------
    Ustar = sqrt( FLXUV )              ! u* [m/s]
    Tstar = -SH / RHOO / CPdry / Ustar ! T* [K]
    Qstar = -LH / RHOO / LHV   / Ustar ! q* [-]
    !Z = ZA - ZDC
    !XXX = 0.4*9.81*Z*TST/TA/UST/UST

    XXX = XXXC
    call cal_psi(XXX,psim,psih)

    XXX2 = (2.0_RP/Z) * XXX
    call cal_psi(XXX2,psim2,psih2)

    XXX10 = (10.0_RP/Z) * XXX
    call cal_psi(XXX10,psim10,psih10)


    FracU10 = ( log(10.0_RP/Z0C ) - psim10 ) / ( log(Z/Z0C ) - psim )
    FracT2  = ( log( 2.0_RP/Z0HC) - psih2  ) / ( log(Z/Z0HC) - psih )

    call BULKFLUX_diagnose_surface( U1, V1,                 & ! [IN]
                                    TA, QA,                 & ! [IN]
                                    RTS, 1.0_RP,            & ! [IN]
                                    Z,                      & ! [IN]
                                    Z0C, Z0HC, 1.0_RP,      & ! [IN]
                                    U10, V10, T2, Q2,       & ! [OUT] ! Q2 is dummy
                                    FracU10, FracT2, 1.0_RP ) ! [IN]

    Q2 = QC

    return
  end subroutine SLC_main

  !-----------------------------------------------------------------------------
  subroutine canopy_wind(ZA, UA, Z0C, ZDC, UC)
    !$acc routine
    implicit none

    real(RP), intent(in)  :: ZA   ! height at 1st atmospheric level [m]
    real(RP), intent(in)  :: UA   ! wind speed at 1st atmospheric level [m/s]
    real(RP), intent(in)  :: Z0C  ! Roughness length above canyon for momentum [m]
    real(RP), intent(in)  :: ZDC  ! Displacement height [m]
    real(RP), intent(out) :: UC   ! wind speed at 1st atmospheric level [m/s]

    real(RP) :: UR,ZC,XLB,BB

    if( ZR + 2.0_RP < ZA ) then
      UR  = UA * log((ZR-ZDC)/Z0C) / log((ZA-ZDC)/Z0C)
      ZC  = 0.7_RP * ZR
      XLB = 0.4_RP * (ZR-ZDC)
      ! BB formulation from Inoue (1963)
      BB  = 0.4_RP * ZR / ( XLB * log((ZR-ZDC)/Z0C) )
      UC  = UR * exp( -BB * (1.0_RP-ZC/ZR) )
    else
      ! PRINT *, 'Warning ZR + 2m  is larger than the 1st WRF level'
      ZC  = ZA / 2.0_RP
      UC  = UA / 2.0_RP
    endif

    UC = max(UC,0.01_RP)

    return
  end subroutine canopy_wind

  !-----------------------------------------------------------------------------
  subroutine cal_beta(BET, BET_CONST, WATER, STRG)
    !$acc routine
    implicit none

    real(RP), intent(out) :: BET       ! evapolation efficiency [-]
    real(RP), intent(in)  :: BET_CONST ! prescribed value
    real(RP), intent(in)  :: WATER     ! rain amount in strage [kg/m2]
    real(RP), intent(in)  :: STRG      ! rain strage [kg/m2]

    if ( BET_CONST > 0.0_RP ) then
       BET = BET_CONST
    else if ( STRG == 0.0_RP ) then ! not consider evapolation from urban
       BET = 0.0_RP
    else
       BET = max( min( WATER / STRG, 1.0_RP ), &
                  1.0E-10_RP ) ! When WATER < STRG/1e10, fix the beta value so that tiny amoumts of water do not coninue to remain
    endif

    return
  end subroutine cal_beta

  !-----------------------------------------------------------------------------
  subroutine cal_psi(zeta,psim,psih)
    !$acc routine
    use scale_const, only: &
       PI     => CONST_PI
    implicit none

    real(RP), intent(inout) :: zeta  ! z/L
    real(RP), intent(out)   :: psim
    real(RP), intent(out)   :: psih
    real(RP)                :: X

    if( zeta >=  1.0_RP ) zeta =  1.0_RP
    if( zeta <= -5.0_RP ) zeta = -5.0_RP

    if( zeta > 0.0_RP ) then
      psim = -5.0_RP * zeta
      psih = -5.0_RP * zeta
    else
      X    = ( 1.0_RP - 16.0_RP * zeta )**0.25_RP
      psim = 2.0_RP * log((1.0_RP+X)/2.0_RP) + log((1.0_RP+X*X)/2.0_RP) - 2.0_RP*atan(X) + PI/2.0_RP
      psih = 2.0_RP * log((1.0_RP+X*X)/2.0_RP)
    end if

    return
  end subroutine cal_psi

  !-----------------------------------------------------------------------------
  !  XXX:   z/L (requires iteration by Newton-Rapson method)
  !  B1:    Stanton number
  !  PSIM:  = PSIX of LSM
  !  PSIH:  = PSIT of LSM
!OCL SERIAL
  subroutine mos(XXX,CH,CD,B1,RIB,Z,Z0,UA,TA,TSF,RHO,i,j)
    !$acc routine
    use scale_const, only: &
       EPS   => CONST_EPS, &
       CPdry => CONST_CPdry ! CPP : heat capacity of dry air [J/K/kg]
    implicit none
integer,intent(in)::i,j
    real(RP), intent(in)    :: B1, Z, Z0, UA, TA, TSF, RHO
    real(RP), intent(out)   :: CD, CH
    real(RP), intent(inout) :: XXX, RIB
    real(RP)                :: XXX0, X, X0, FAIH, DPSIM, DPSIH
    real(RP)                :: F, DF, XXXP, US, TS, AL, XKB, DD, PSIM, PSIH
    integer                 :: NEWT
    integer, parameter      :: NEWT_END = 10

    real(RP)                :: lnZ
    real(RP)                :: sqX, sqX0

    lnZ = log( (Z+Z0)/Z0 )

    if( RIB < 0.0_RP ) then

       RIB = max( RIB, -15.0_RP )

       do NEWT = 1, NEWT_END

          XXX = min( XXX, -1.0E-3_RP )

          XXX0  = XXX * Z0/(Z+Z0)

          sqX   = sqrt( 1.0_RP - 16.0_RP * XXX  )
          sqX0  = sqrt( 1.0_RP - 16.0_RP * XXX0 )

          X     = sqrt( sqX )
          X0    = sqrt( sqX0 )

          PSIM  = lnZ &
                - log( (X+1.0_RP)**2 * (sqX  + 1.0_RP) ) &
                + 2.0_RP * atan(X) &
                + log( (X+1.0_RP)**2 * (sqX0 + 1.0_RP) ) &
                - 2.0_RP * atan(X0)
          FAIH  = 1.0_RP / sqX
          PSIH  = lnZ + 0.4_RP*B1 &
                - 2.0_RP * log( sqX  + 1.0_RP ) &
                + 2.0_RP * log( sqX0 + 1.0_RP )

          DPSIM = 1.0_RP / ( X  * XXX ) &
                - 1.0_RP / ( X0 * XXX )
          DPSIH = 1.0_RP / ( sqX  * XXX ) &
                - 1.0_RP / ( sqX0 * XXX )

          F     = RIB * PSIM**2 / PSIH - XXX

          DF    = RIB * ( 2.0_RP*DPSIM*PSIM*PSIH - DPSIH*PSIM**2 ) &
                / PSIH**2 - 1.0_RP

          XXXP  = XXX
          XXX   = XXXP - F / DF

          XXX = max( XXX, -10.0_RP )

       end do

    else if( RIB >= 0.142857_RP ) then

       XXX  = 0.714_RP
       PSIM = lnZ + 7.0_RP * XXX
       PSIH = PSIM + 0.4_RP * B1

    else

       AL   = lnZ
       XKB  = 0.4_RP * B1
       DD   = -4.0_RP * RIB * 7.0_RP * XKB * AL + (AL+XKB)**2
       DD = max( DD, 0.0_RP )

       XXX  = ( AL + XKB - 2.0_RP*RIB*7.0_RP*AL - sqrt(DD) ) / ( 2.0_RP * ( RIB*7.0_RP**2 - 7.0_RP ) )
       XXX = min( XXX, 0.714_RP )
       PSIM = lnZ + 7.0_RP * XXX
       PSIH = PSIM + 0.4_RP * B1

    endif

    US = 0.4_RP * UA / PSIM             ! u*
    if( US <= 0.01_RP ) US = 0.01_RP
    !TS = 0.4_RP * (TA-TSF) / PSIH       ! T*

    CD    = US * US / (UA+EPS)**2         ! CD
    CH    = 0.4_RP * US / PSIH / (UA+EPS) ! CH
    !ALPHA = RHO * CPdry * 0.4_RP * US / PSIH  ! RHO*CP*CH*U

    return
  end subroutine mos

  !-------------------------------------------------------------------
!OCL SERIAL
  subroutine multi_layer( &
       KM,BOUND, &
       A, B, C, D, P, Q, &
       G0,CAP,AKS,TSL,DZ,DELT,TSLEND)
  !
  !  calculate temperature in roof/building/road
  !  multi-layer heat equation model
  !  Solving Heat Equation by Tri Diagonal Matrix Algorithm
  !-------------------------------------------------------------------
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: KM
    integer,  intent(in)    :: BOUND
    real(RP), intent(in)    :: G0
    real(RP), intent(in)    :: CAP
    real(RP), intent(in)    :: AKS
    real(RP), intent(inout) :: TSL(KM)
    real(RP), intent(in)    :: DZ(KM)
    real(RP), intent(in)    :: DELT      ! Tim setep [ s ]
    real(RP), intent(in)    :: TSLEND
    real(RP), intent(out)   :: A(KM), B(KM), C(KM), D(KM), P(KM), Q(KM)

    real(RP)                :: DZEND
    integer                 :: K

    DZEND = DZ(KM)

    A(1)  = 0.0_RP

    B(1)  = CAP * DZ(1) / DELT &
          + 2.0_RP * AKS / (DZ(1)+DZ(2))
    C(1)  = -2.0_RP * AKS / (DZ(1)+DZ(2))
    D(1)  = CAP * DZ(1) / DELT * TSL(1) + G0

    do K = 2, KM-1
      A(K) = -2.0_RP * AKS / (DZ(K-1)+DZ(K))
      B(K) = CAP * DZ(K) / DELT + 2.0_RP * AKS / (DZ(K-1)+DZ(K)) + 2.0_RP * AKS / (DZ(K)+DZ(K+1))
      C(K) = -2.0_RP * AKS / (DZ(K)+DZ(K+1))
      D(K) = CAP * DZ(K) / DELT * TSL(K)
    end do

    if( BOUND == 1 ) then ! Flux=0
      A(KM) = -2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      B(KM) = CAP * DZ(KM) / DELT + 2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      C(KM) = 0.0_RP
      D(KM) = CAP * DZ(KM) / DELT * TSL(KM)
    else ! T=constant
      A(KM) = -2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      B(KM) = CAP * DZ(KM) / DELT + 2.0_RP * AKS / (DZ(KM-1)+DZ(KM)) + 2.0_RP * AKS / (DZ(KM)+DZEND)
      C(KM) = 0.0_RP
      D(KM) = CAP * DZ(KM) / DELT * TSL(KM) + 2.0_RP * AKS * TSLEND / (DZ(KM)+DZEND)
    end if

    P(1) = -C(1) / B(1)
    Q(1) =  D(1) / B(1)

    !$acc loop seq
    do K = 2, KM
      P(K) = -C(K) / ( A(K) * P(K-1) + B(K) )
      Q(K) = ( -A(K) * Q(K-1) + D(K) ) / ( A(K) * P(K-1) + B(K) )
    end do

    TSL(KM) = Q(KM)

    !$acc loop seq
    do K = KM-1, 1, -1
      TSL(K) = P(K) * TSL(K+1) + Q(K)
    end do

    return
  end subroutine multi_layer

!!$  !-------------------------------------------------------------------
!!$!OCL SERIAL
!!$  subroutine multi_layer2(KM,BOUND,G0,CAP,AKS,TSL,DZ,DELT,TSLEND,CAP1,AKS1)
!!$  !
!!$  !  calculate temperature in roof/building/road
!!$  !  multi-layer heat equation model
!!$  !  Solving Heat Equation by Tri Diagonal Matrix Algorithm
!!$  !-------------------------------------------------------------------
!!$
!!$    implicit none
!!$
!!$    real(RP), intent(in)    :: G0
!!$    real(RP), intent(in)    :: CAP
!!$    real(RP), intent(in)    :: AKS
!!$    real(RP), intent(in)    :: CAP1      ! for 1st layer
!!$    real(RP), intent(in)    :: AKS1      ! for 1st layer
!!$    real(DP), intent(in)    :: DELT      ! Time step [ s ]
!!$    real(RP), intent(in)    :: TSLEND
!!$    integer,  intent(in)    :: KM
!!$    integer,  intent(in)    :: BOUND
!!$    real(RP), intent(in)    :: DZ(KM)
!!$    real(RP), intent(inout) :: TSL(KM)
!!$    real(RP)                :: A(KM), B(KM), C(KM), D(KM), X(KM), P(KM), Q(KM)
!!$    real(RP)                :: DZEND
!!$    integer                 :: K
!!$
!!$    DZEND = DZ(KM)
!!$
!!$    A(1)  = 0.0_RP
!!$
!!$    B(1)  = CAP1 * DZ(1) / DELT &
!!$          + 2.0_RP * AKS1 / (DZ(1)+DZ(2))
!!$    C(1)  = -2.0_RP * AKS1 / (DZ(1)+DZ(2))
!!$    D(1)  = CAP1 * DZ(1) / DELT * TSL(1) + G0
!!$
!!$    do K = 2, KM-1
!!$      A(K) = -2.0_RP * AKS / (DZ(K-1)+DZ(K))
!!$      B(K) = CAP * DZ(K) / DELT + 2.0_RP * AKS / (DZ(K-1)+DZ(K)) + 2.0_RP * AKS / (DZ(K)+DZ(K+1))
!!$      C(K) = -2.0_RP * AKS / (DZ(K)+DZ(K+1))
!!$      D(K) = CAP * DZ(K) / DELT * TSL(K)
!!$    end do
!!$
!!$    if( BOUND == 1 ) then ! Flux=0
!!$      A(KM) = -2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
!!$      B(KM) = CAP * DZ(KM) / DELT + 2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
!!$      C(KM) = 0.0_RP
!!$      D(KM) = CAP * DZ(KM) / DELT * TSL(KM)
!!$    else ! T=constant
!!$      A(KM) = -2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
!!$      B(KM) = CAP * DZ(KM) / DELT + 2.0_RP * AKS / (DZ(KM-1)+DZ(KM)) + 2.0_RP * AKS / (DZ(KM)+DZEND)
!!$      C(KM) = 0.0_RP
!!$      D(KM) = CAP * DZ(KM) / DELT * TSL(KM) + 2.0_RP * AKS * TSLEND / (DZ(KM)+DZEND)
!!$    end if
!!$
!!$    P(1) = -C(1) / B(1)
!!$    Q(1) =  D(1) / B(1)
!!$
!!$    !$acc loop seq
!!$    do K = 2, KM
!!$      P(K) = -C(K) / ( A(K) * P(K-1) + B(K) )
!!$      Q(K) = ( -A(K) * Q(K-1) + D(K) ) / ( A(K) * P(K-1) + B(K) )
!!$    end do
!!$
!!$    X(KM) = Q(KM)
!!$
!!$    !$acc loop seq
!!$    do K = KM-1, 1, -1
!!$      X(K) = P(K) * X(K+1) + Q(K)
!!$    end do
!!$
!!$    do K = 1, KM
!!$      TSL(K) = X(K)
!!$    enddo
!!$
!!$    return
!!$  end subroutine multi_layer2

  !-----------------------------------------------------------------------------
  !> set urban parameters
  subroutine urban_param_setup
    implicit none

    real(RP) :: DHGT,THGT,VFWS,VFGS
    integer  :: k

    ! initialize
    R        = 0.0_RP
    RW       = 0.0_RP
    HGT      = 0.0_RP
    Z0HR     = 0.0_RP
    Z0HB     = 0.0_RP
    Z0HG     = 0.0_RP
    ZDC_TBL  = 0.0_RP
    Z0C_TBL  = 0.0_RP
    Z0HC_TBL = 0.0_RP
    SVF      = 0.0_RP

    ! set up other urban parameters
    Z0HR     = 0.1_RP * Z0R
    Z0HB     = 0.1_RP * Z0B
    Z0HG     = 0.1_RP * Z0G
    ZDC_TBL  = ZR * 0.3_RP
    Z0C_TBL  = ZR * 0.15_RP
    Z0HC_TBL = 0.1_RP * Z0C_TBL

    ! HGT:  Normalized height
    HGT  = ZR / ( ROAD_WIDTH + ROOF_WIDTH )

    ! R:  Normalized Roof Width (a.k.a. "building coverage ratio")
    R    = ROOF_WIDTH / ( ROAD_WIDTH + ROOF_WIDTH )
    RW   = 1.0_RP - R

    ! Calculate Sky View Factor:
    DHGT = HGT / 100.0_RP
    THGT = 0.0_RP
    VFWS = 0.0_RP
    THGT = HGT - DHGT / 2.0_RP
    do k = 1, 99
      THGT  = THGT - DHGT
      VFWS = VFWS + 0.25_RP * ( 1.0_RP - THGT / sqrt( THGT**2 + RW**2 ) )
    end do

    VFWS = VFWS / 99.0_RP
    VFWS = VFWS * 2.0_RP
    VFGS = 1.0_RP - 2.0_RP * VFWS * HGT / RW
    SVF  = VFGS

    return
  end subroutine urban_param_setup

  !-----------------------------------------------------------------------------
  !> read urban data from paraneter table
  subroutine read_urban_param_table( INFILENAME )
    use scale_prc, only: &
        PRC_abort
    implicit none

    character(*), intent(in) :: INFILENAME

    real(RP) :: BETR, BETB, BETG

    namelist / PARAM_URBAN_DATA / &
       ZR,         &
       roof_width, &
       road_width, &
       SIGMA_ZED,  &
       AH_TBL,     &
       AHL_TBL,    &
       BETR,       &
       BETB,       &
       BETG,       &
       STRGR,      &
       STRGB,      &
       STRGG,      &
       CAPR,       &
       CAPB,       &
       CAPG,       &
       AKSR,       &
       AKSB,       &
       AKSG,       &
       ALBR,       &
       ALBB,       &
       ALBG,       &
       EPSR,       &
       EPSB,       &
       EPSG,       &
       Z0R,        &
       Z0B,        &
       Z0G,        &
       TRLEND,     &
       TBLEND,     &
       TGLEND

    character(len=H_LONG) :: fname
    integer :: fid
    integer :: ierr
    !---------------------------------------------------------------------------

    fid = IO_get_available_fid()
    call IO_get_fname(fname, INFILENAME)
    open( fid,                  &
          file   = fname,       &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )
    if ( ierr /= 0 ) then
       LOG_NEWLINE
       LOG_ERROR("URBAN_DYN_kusaka01_setup",*) 'Failed to open a parameter file : ', trim(fname)
       call PRC_abort
    end if

    LOG_NEWLINE
    LOG_INFO("URBAN_DYN_kusaka01_setup",*) 'read_urban_param_table: Read urban parameters from file'

    BETR = -1.0_RP
    BETB = -1.0_RP
    BETG = -1.0_RP

    !--- read namelist
    rewind(fid)
    read  (fid,nml=PARAM_URBAN_DATA,iostat=ierr)
    if ( ierr < 0 ) then !--- no data
       LOG_INFO("read_urban_param_table",*)  'Not found namelist of PARAM_URBAN_DATA. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("read_urban_param_table",*) 'Not appropriate names in PARAM_URBAN_DATA of ', trim(INFILENAME)
       call PRC_abort
    endif
    LOG_NML(PARAM_URBAN_DATA)

    close( fid )

    if ( ZR <= 0.0_RP ) then
       LOG_ERROR("read_urban_param_table",*) 'ZR is not appropriate value; ZR must be larger than 0. ZR=', ZR
       call PRC_abort
    endif

    BETR_CONST = BETR
    BETB_CONST = BETB
    BETG_CONST = BETG

    return
  end subroutine read_urban_param_table

  !-----------------------------------------------------------------------------
  !> read 2D gridded urban data
  subroutine read_urban_gridded_data_2D(  &
       UIA, UJA,                          &
       INFILENAME,                        &
       VARNAME,                           &
       udata                              )
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush, &
       FILE_CARTESC_check_coordinates, &
       FILE_CARTESC_close
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer         , intent(in)  :: UIA, UJA
    character(len=*), intent(in)  :: INFILENAME
    character(len=*), intent(in)  :: VARNAME
    real(RP),         intent(out) :: udata(UIA,UJA) !< value of the variable

    integer :: fid
    !---------------------------------------------------------------------------
    LOG_NEWLINE
    LOG_INFO("URBAN_DYN_kusaka01_setup",*) 'read_urban_gridded_data ',trim(VARNAME)

    fid = IO_get_available_fid()
    call FILE_CARTESC_open( INFILENAME, fid )
    call FILE_CARTESC_read( fid, VARNAME, 'XY', udata(:,:) )

    call FILE_CARTESC_flush( fid )
    call FILE_CARTESC_check_coordinates( fid )
    call FILE_CARTESC_close( fid )

    return
  end subroutine read_urban_gridded_data_2D

  !-----------------------------------------------------------------------------
  !> read 3D gridded urban data
  subroutine read_urban_gridded_data_3D(  &
       UIA, UJA,                          &
       INFILENAME,                        &
       VARNAME,                           &
       udata                              )
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush, &
       FILE_CARTESC_check_coordinates, &
       FILE_CARTESC_close
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer         , intent(in)  :: UIA, UJA
    character(len=*), intent(in)  :: INFILENAME
    character(len=*), intent(in)  :: VARNAME
    real(RP),         intent(out) :: udata(UIA,UJA,24) !< value of the variable

    integer :: fid, k
    !---------------------------------------------------------------------------
    LOG_NEWLINE
    LOG_INFO("URBAN_DYN_kusaka01_setup",*) 'read_urban_gridded_data ',trim(VARNAME)

    fid = IO_get_available_fid()
    call FILE_CARTESC_open( INFILENAME, fid )
    do k = 1, 24
    call FILE_CARTESC_read( fid, VARNAME, 'XY', udata(:,:,k), k )
    enddo

    call FILE_CARTESC_flush( fid )
    call FILE_CARTESC_check_coordinates( fid )
    call FILE_CARTESC_close( fid )

    return
  end subroutine read_urban_gridded_data_3D

  !-----------------------------------------------------------------------------
  subroutine put_history( &
       UIA, UJA, &
       SHR, SHB, SHG, &
       LHR, LHB, LHG, &
       GHR, GHB, GHG, &
       RNR, RNB, RNG, &
       RNgrd          )
    use scale_file_history, only: &
       FILE_HISTORY_put
    integer, intent(in) :: UIA, UJA
    real(RP), intent(in) :: SHR(UIA,UJA), SHB(UIA,UJA), SHG(UIA,UJA)
    real(RP), intent(in) :: LHR(UIA,UJA), LHB(UIA,UJA), LHG(UIA,UJA)
    real(RP), intent(in) :: GHR(UIA,UJA), GHB(UIA,UJA), GHG(UIA,UJA)
    real(RP), intent(in) :: RNR(UIA,UJA), RNB(UIA,UJA), RNG(UIA,UJA)
    real(RP), intent(in) :: RNgrd(UIA,UJA)

    call FILE_HISTORY_put( I_SHR,   SHR  (:,:) )
    call FILE_HISTORY_put( I_SHB,   SHB  (:,:) )
    call FILE_HISTORY_put( I_SHG,   SHG  (:,:) )
    call FILE_HISTORY_put( I_LHR,   LHR  (:,:) )
    call FILE_HISTORY_put( I_LHB,   LHB  (:,:) )
    call FILE_HISTORY_put( I_LHG,   LHG  (:,:) )
    call FILE_HISTORY_put( I_GHR,   GHR  (:,:) )
    call FILE_HISTORY_put( I_GHB,   GHB  (:,:) )
    call FILE_HISTORY_put( I_GHG,   GHG  (:,:) )
    call FILE_HISTORY_put( I_RNR,   RNR  (:,:) )
    call FILE_HISTORY_put( I_RNB,   RNB  (:,:) )
    call FILE_HISTORY_put( I_RNG,   RNG  (:,:) )
    call FILE_HISTORY_put( I_RNgrd, RNgrd(:,:) )

    return
  end subroutine put_history

end module scale_urban_dyn_kusaka01
