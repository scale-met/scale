!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!      Put atmospheric data for urban test
!!      Test is based on Kusaka et al. (2000,BLM)
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
       I_LW  => CONST_I_LW,  &
       I_SW  => CONST_I_SW,  &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, &
       PRE00 => CONST_PRE00, &
       LHV   => CONST_LHV,   &    ! ELL : latent heat of vaporization [J/kg]
       Rdry  => CONST_Rdry,  &    ! specific gas constant (dry air)
       CPdry => CONST_CPdry       ! specific heat (dry air,constant pressure) [J/kg/K]
    use scale_grid_real, only: &
       REAL_lon
    use mod_cpl_vars, only: &
       TMPA  => URB_ATM_TEMP,        &
       PRSA  => URB_ATM_PRES,        &
       WA    => URB_ATM_W,           &
       UA    => URB_ATM_U,           &
       VA    => URB_ATM_V,           &
       RHOA  => URB_ATM_DENS,        &
       QVA   => URB_ATM_QV,          &
       PBL   => URB_ATM_PBL,         &
       PRSS  => URB_ATM_SFC_PRES,    &
       RWD   => URB_ATM_SFLX_rad_dn, &
       RAIN  => URB_ATM_SFLX_rain,   &
       SNOW  => URB_ATM_SFLX_snow
    use scale_time, only:   &
       NOWSEC => TIME_NOWSEC,      & !< subday part  of current time [sec]
       dt_URB => TIME_DTSEC_URBAN    !< time interval of urban step  [sec]
    use mod_admin_time, only: &
       TIME_DOURBAN_step             !< execute urban component in this step?
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP) :: WORK(IA,JA)
    real(RP) :: LON, LAT
    real(RP) :: dsec
    integer  :: tloc

    real(RP) :: LWD(IA,JA)
    real(RP) :: SWD(IA,JA)

    real(RP), parameter :: SRATIO = 0.75_RP ! ratio between direct/total solar [-]
    real(RP)            :: SWtot

    real(RP) :: SW  (0:24)
    real(RP) :: PT  (0:24)
    real(RP) :: Wind(0:24)

    data SW /     0.0_RP, & ! 00LST
                  0.0_RP, & ! 01LST
                  0.0_RP, & ! 02LST
                  0.0_RP, & ! 03LST
                  0.0_RP, & ! 04LST
                  0.0_RP, & ! 05LST
                 50.0_RP, & ! 06LST
                240.0_RP, & ! 07LST
                420.0_RP, & ! 08LST
                600.0_RP, & ! 09LST
                690.0_RP, & ! 10LST
                765.0_RP, & ! 11LST
                800.0_RP, & ! 12LST
                765.0_RP, & ! 13LST
                690.0_RP, & ! 14LST
                600.0_RP, & ! 15LST
                420.0_RP, & ! 16LST
                240.0_RP, & ! 17LST
                 50.0_RP, & ! 18LST
                  0.0_RP, & ! 19LST
                  0.0_RP, & ! 20LST
                  0.0_RP, & ! 21LST
                  0.0_RP, & ! 22LST
                  0.0_RP, & ! 23LST
                  0.0_RP  / ! 24LST

    data PT /   28.00_RP, & ! 00LST
                27.80_RP, & ! 01LST
                27.65_RP, & ! 02LST
                27.50_RP, & ! 03LST
                27.35_RP, & ! 04LST
                27.20_RP, & ! 05LST
                27.10_RP, & ! 06LST
                27.50_RP, & ! 07LST
                27.85_RP, & ! 08LST
                28.25_RP, & ! 09LST
                28.80_RP, & ! 10LST
                29.40_RP, & ! 11LST
                30.00_RP, & ! 12LST
                30.20_RP, & ! 13LST
                30.40_RP, & ! 14LST
                30.60_RP, & ! 15LST
                30.35_RP, & ! 16LST
                30.10_RP, & ! 17LST
                29.85_RP, & ! 18LST
                29.55_RP, & ! 19LST
                29.15_RP, & ! 20LST
                28.75_RP, & ! 21LST
                28.50_RP, & ! 22LST
                28.25_RP, & ! 23LST
                28.00_RP  / ! 24LST

    data Wind /  2.75_RP, & ! 00LST
                 2.75_RP, & ! 01LST
                 2.75_RP, & ! 02LST
                 2.75_RP, & ! 03LST
                 2.75_RP, & ! 04LST
                 2.75_RP, & ! 05LST
                 2.80_RP, & ! 06LST
                 3.00_RP, & ! 07LST
                 3.25_RP, & ! 08LST
                 3.50_RP, & ! 09LST
                 3.65_RP, & ! 10LST
                 3.65_RP, & ! 11LST
                 3.50_RP, & ! 12LST
                 3.40_RP, & ! 13LST
                 3.27_RP, & ! 14LST
                 3.15_RP, & ! 15LST
                 3.05_RP, & ! 16LST
                 2.95_RP, & ! 17LST
                 2.85_RP, & ! 18LST
                 2.80_RP, & ! 19LST
                 2.75_RP, & ! 20LST
                 2.70_RP, & ! 21LST
                 2.72_RP, & ! 22LST
                 2.75_RP, & ! 23LST
                 2.75_RP  / ! 24LST

    real(RP) :: THETA
    real(RP) :: RovCP
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    RovCP = Rdry / CPdry

    if ( USER_do ) then

       VA  (:,:)        =      0.0_RP
       WA  (:,:)        =      0.0_RP
       RHOA(:,:)        =     1.13_RP
       PBL (:,:)        =    100.0_RP
       RWD (:,:,I_LW,1) =      0.0_RP ! direct
       RWD (:,:,I_LW,2) =    400.0_RP ! diffuse
       PRSA(:,:)        = 100000.0_RP
       PRSS(:,:)        = 100120.0_RP
       QVA (:,:)        =    0.015_RP
       RAIN(:,:)        =      0.0_RP
       SNOW(:,:)        =      0.0_RP

       do j = 1, JA
       do i = 1, IA

          LON = REAL_lon(i,j) / D2R

          tloc = mod(int(NOWSEC/3600.0_RP)+int(LON/15.0_RP),24)

          dsec = mod(NOWSEC,3600.0_RP) / 3600.0_RP

          SWtot = ( ( 1.0_RP-dsec ) * SW(tloc  ) &
                  + (        dsec ) * SW(tloc+1) )

          RWD (i,j,I_SW,1) = (        SRATIO ) * SWtot ! direct
          RWD (i,j,I_SW,2) = ( 1.0_RP-SRATIO ) * SWtot ! diffuse

          THETA     = ( ( 1.0_RP-dsec ) * PT(tloc  ) &
                      + (        dsec ) * PT(tloc+1) ) + TEM00
          TMPA(i,j) = THETA * ( PRSA(i,j) / PRE00 )**RovCP   ! air temp, but now PRSA = 100000Pa

          UA  (i,j) = ( ( 1.0_RP-dsec ) * Wind(tloc  ) &
                      + (        dsec ) * Wind(tloc+1) )

          LWD (i,j) = RWD(i,j,I_LW,1) + RWD(i,j,I_LW,2)
          SWD (i,j) = RWD(i,j,I_SW,1) + RWD(i,j,I_SW,2)

       enddo
       enddo

       call HIST_in( WORK (:,:), 'PT_urb',   'Potential air temperature',    'K'     )
       call HIST_in( QVA  (:,:), 'QA_urb',   'Specific humidity',            'kg/kg' )
       call HIST_in( UA   (:,:), 'UA_urb',   'Wind speed',                   'm/s'   )
       call HIST_in( SWD  (:,:), 'SWD_urb',  'Downward shortwave radiation', 'W/m2'  )
       call HIST_in( LWD  (:,:), 'LWD_urb',  'Downward longwave radiation',  'W/m2'  )

       WORK(:,:) = ( RAIN(:,:) + SNOW(:,:) ) * dt_URB
       call HIST_in( WORK (:,:), 'RAIN_urb', 'Precipitation',                'kg/m2' )

    endif

    return
  end subroutine USER_step

end module mod_user
