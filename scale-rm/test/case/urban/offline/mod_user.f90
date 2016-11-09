!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!      Put atmospheric data for urban test
!!      Test is based on Aoyagi et al. (2011,JAMC)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

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
  logical, private :: USER_do   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_atmos_admin, only: &
       ATMOS_do,        &
       ATMOS_sw_dyn,    &
       ATMOS_sw_phy_mp, &
       ATMOS_sw_phy_ae, &
       ATMOS_sw_phy_ch, &
       ATMOS_sw_phy_rd, &
       ATMOS_sw_phy_sf, &
       ATMOS_sw_phy_tb, &
       ATMOS_sw_phy_cp
    implicit none

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_USER)

    if( IO_L ) write(IO_FID_LOG,*)

    ! atmosphric model set to off
    ATMOS_do        = .false.
    ATMOS_sw_dyn    = .false.
    ATMOS_sw_phy_mp = .false.
    ATMOS_sw_phy_ae = .false.
    ATMOS_sw_phy_ch = .false.
    ATMOS_sw_phy_rd = .false.
    ATMOS_sw_phy_sf = .false.
    ATMOS_sw_phy_tb = .false.
    ATMOS_sw_phy_cp = .false.

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    ! replace urban variables and input to history buffer
    call USER_step

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_const, only: &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, &
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &    ! specific gas constant (dry air)
       CPdry => CONST_CPdry       ! specific heat (dry air,constant pressure) [J/kg/K]
    use scale_grid_real, only: &
       REAL_lon
    use mod_cpl_vars, only: &
       TMPA  => URB_ATM_TEMP,      &  ! air temperature
       PRSA  => URB_ATM_PRES,      &
       WA    => URB_ATM_W,         &
       UA    => URB_ATM_U,         &
       VA    => URB_ATM_V,         &
       RHOA  => URB_ATM_DENS,      &
       QVA   => URB_ATM_QV,        &
       PBL   => URB_ATM_PBL,       &
       PRSS  => URB_ATM_SFC_PRES,  &
       LWD   => URB_ATM_SFLX_LW,   &
       SWD   => URB_ATM_SFLX_SW,   &
       PREC  => URB_ATM_SFLX_prec
    use mod_urban_vars, only: &
       URBAN_TR,         &
       URBAN_TB,         &
       URBAN_TG,         &
       URBAN_TC,         &
       URBAN_QC,         &
       URBAN_UC,         &
       URBAN_TRL,        &
       URBAN_TBL,        &
       URBAN_TGL,        &
       URBAN_ROFF
    use scale_time, only:   &
       NOWSEC => TIME_NOWSEC,      & !< absolute sec
       dt_URB => TIME_DTSEC_URBAN, & !< time interval of urban step [sec]
       TIME_DOURBAN_step
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP) :: WORK(IA,JA)
    real(RP) :: LON, LAT
    real(RP) :: dsec
    integer  :: tloc

    real(RP) :: SW
    real(RP) :: PT
    real(RP) :: Wind
    real(RP) :: Rain
    real(RP) :: Qvapor ! mixing ratio [kg/kg]
    real(RP) :: ES

    real(RP) :: THETA

    real(RP) :: SX,SB,SG
    real(RP) :: VFGS,VFWS,VFGW,VFWG

    real(RP) :: RovCP
    integer :: k, i, j
    !---------------------------------------------------------------------------

    RovCP = Rdry / CPdry

    if ( USER_do ) then

       UA  (:,:) = 18.95165494678601_RP
       VA  (:,:)  = 0.0_RP
       WA  (:,:)  = 0.0_RP
       RHOA(:,:)  = 1.193221659609323_RP
       PBL (:,:)  = 100.0_RP
       LWD (:,:)  = 434.6034964144717_RP
       PRSA(:,:)  = 100000.0_RP
       PRSS(:,:)  = 102400.6750905938_RP
       QVA (:,:)  = 1.612903266525567E-02_RP
       PREC(:,:)  = 5.0_RP / 3600.0_RP

       URBAN_ROFF(:,:) = 0.0_RP
       !URBAN_TR(:,:) = 300.0_RP
       URBAN_TB(:,:) = 298.4421947975652_RP
       URBAN_TG(:,:) = 298.6575894549702_RP
       URBAN_TC(:,:) = 297.7168673290382_RP
       URBAN_QC(:,:) = 1.830501916072291E-02_RP
       URBAN_UC(:,:) = 7.978626675770874_RP
       !
       URBAN_TBL(1,:,:) = 298.4421947975652_RP
       URBAN_TBL(2,:,:) = 298.9831365435758_RP
       URBAN_TBL(3,:,:) = 299.6502356623731_RP
       URBAN_TBL(4,:,:) = 299.6819173427097_RP
       URBAN_TBL(5,:,:) = 299.0036771970256_RP
       URBAN_TGL(1,:,:) = 298.6575894549702_RP
       URBAN_TGL(2,:,:) = 299.2942430001926_RP
       URBAN_TGL(3,:,:) = 300.0760515131021_RP
       URBAN_TGL(4,:,:) = 300.0731793271447_RP
       URBAN_TGL(5,:,:) = 299.1892611738443_RP

       SB   = 1.429919681745362_RP
       SG   = 2.090724511380252_RP
       VFGS = 0.5335010498145294_RP
       VFWS = 0.3498742126391030_RP
       VFGW = 0.4664989501854706_RP
       VFWG = 0.3498742126391030_RP

       SX = SG / ( VFGS * 0.8_RP +  VFWS * 0.8_RP * 0.2_RP / 0.8_RP * VFGW * 0.8_RP  )
       print *,"SX1",SX
       SX = SB / ( VFWS * 0.8_RP +  VFGS * 0.8_RP * 0.2_RP / 0.8_RP * VFWG * 0.8_RP  )
       print *,"SX2",SX
       SWD (:,:) = SX

       do j = 1, JA
       do i = 1, IA

       !   LON = REAL_lon(i,j) / D2R
       !
       !   tloc = mod(int(NOWSEC/3600.0_RP)+int(LON/15.0_RP),24)
       !
       !   dsec = mod(NOWSEC,3600.0_RP) / 3600.0_RP
       !
           THETA = 293.7453140572144_RP
           WORK(i,j) = THETA
           TMPA(i,j) = THETA / ( PRE00 / PRSA(i,j) )**RovCP

       enddo
       enddo

       call HIST_in( WORK (:,:), 'PT_urb',   'Potential air temperature',    'K'     )
       call HIST_in( QVA  (:,:), 'QA_urb',   'Specific humidity',            'kg/kg' )
       call HIST_in( UA   (:,:), 'UA_urb',   'Wind speed',                   'm/s'   )
       call HIST_in( SWD  (:,:), 'SWD_urb',  'Downward shortwave radiation', 'W/m2'  )
       call HIST_in( LWD  (:,:), 'LWD_urb',  'Downward longwave radiation',  'W/m2'  )

       WORK(:,:) = PREC(:,:) * dt_URB
       call HIST_in( WORK (:,:), 'RAIN_urb', 'Precipitation',                'kg/m2' )

    endif

    return
  end subroutine USER_step

end module mod_user
