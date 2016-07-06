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
       PRE00 => CONST_PRE00, &
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
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP) :: WORK(IA,JA)
    real(RP) :: LON
    real(RP) :: dsec
    integer  :: tloc

    real(RP) :: LWD(IA,JA)
    real(RP) :: SWD(IA,JA)

    real(RP), parameter :: SRATIO = 0.75_RP ! ratio between direct/total solar [-]
    real(RP)            :: SWtot

    real(RP) :: SW    (0:24)
    real(RP) :: PT    (0:24)
    real(RP) :: Wind  (0:24)
    real(RP) :: Prcp  (0:24)
    real(RP) :: Qvapor(0:24) ! mixing ratio [kg/kg]

    data SW /       0.0_RP, & ! 00LST
                    0.0_RP, & ! 01LST
                    0.0_RP, & ! 02LST
                    0.0_RP, & ! 03LST
                    0.0_RP, & ! 04LST
                   11.0_RP, & ! 05LST
                   26.5_RP, & ! 06LST
                  161.5_RP, & ! 07LST
                  434.2_RP, & ! 08LST
                  694.6_RP, & ! 09LST
                  877.5_RP, & ! 10LST
                  807.4_RP, & ! 11LST
                  907.4_RP, & ! 12LST
                  847.2_RP, & ! 13LST
                  725.2_RP, & ! 14LST
                  559.0_RP, & ! 15LST
                  357.9_RP, & ! 16LST
                  158.6_RP, & ! 17LST
                   14.7_RP, & ! 18LST
                    0.0_RP, & ! 19LST
                    0.0_RP, & ! 20LST
                    0.0_RP, & ! 21LST
                    0.0_RP, & ! 22LST
                    0.0_RP, & ! 23LST
                    0.0_RP  / ! 24LST

    data PT /     296.6_RP, & ! 00LST
                  296.6_RP, & ! 01LST
                  296.7_RP, & ! 02LST
                  296.5_RP, & ! 03LST
                  296.4_RP, & ! 04LST
                  296.3_RP, & ! 05LST
                  295.9_RP, & ! 06LST
                  296.1_RP, & ! 07LST
                  296.7_RP, & ! 08LST
                  298.0_RP, & ! 09LST
                  299.2_RP, & ! 10LST
                  300.0_RP, & ! 11LST
                  300.6_RP, & ! 12LST
                  301.7_RP, & ! 13LST
                  302.2_RP, & ! 14LST
                  302.4_RP, & ! 15LST
                  302.0_RP, & ! 16LST
                  301.1_RP, & ! 17LST
                  300.5_RP, & ! 18LST
                  299.8_RP, & ! 19LST
                  299.1_RP, & ! 20LST
                  298.9_RP, & ! 21LST
                  298.6_RP, & ! 22LST
                  298.3_RP, & ! 23LST
                  296.6_RP  / ! 24LST

    data Wind /     3.4_RP, & ! 00LST
                    3.2_RP, & ! 01LST
                    3.2_RP, & ! 02LST
                    3.2_RP, & ! 03LST
                    3.2_RP, & ! 04LST
                    3.0_RP, & ! 05LST
                    3.8_RP, & ! 06LST
                    3.6_RP, & ! 07LST
                    3.5_RP, & ! 08LST
                    3.6_RP, & ! 09LST
                    4.0_RP, & ! 10LST
                    3.8_RP, & ! 11LST
                    4.0_RP, & ! 12LST
                    3.9_RP, & ! 13LST
                    4.0_RP, & ! 14LST
                    4.1_RP, & ! 15LST
                    3.7_RP, & ! 16LST
                    3.6_RP, & ! 17LST
                    4.2_RP, & ! 18LST
                    3.9_RP, & ! 19LST
                    3.7_RP, & ! 20LST
                    3.4_RP, & ! 21LST
                    3.8_RP, & ! 22LST
                    4.0_RP, & ! 23LST
                    3.4_RP  / ! 24LST

    data Qvapor / 0.018_RP, & ! 00LST
                  0.018_RP, & ! 01LST
                  0.018_RP, & ! 02LST
                  0.018_RP, & ! 03LST
                  0.018_RP, & ! 04LST
                  0.018_RP, & ! 05LST
                  0.017_RP, & ! 06LST
                  0.017_RP, & ! 07LST
                  0.017_RP, & ! 08LST
                  0.018_RP, & ! 09LST
                  0.018_RP, & ! 10LST
                  0.018_RP, & ! 11LST
                  0.018_RP, & ! 12LST
                  0.017_RP, & ! 13LST
                  0.017_RP, & ! 14LST
                  0.016_RP, & ! 15LST
                  0.017_RP, & ! 16LST
                  0.018_RP, & ! 17LST
                  0.018_RP, & ! 18LST
                  0.017_RP, & ! 19LST
                  0.018_RP, & ! 20LST
                  0.018_RP, & ! 21LST
                  0.018_RP, & ! 22LST
                  0.018_RP, & ! 23LST
                  0.018_RP  / ! 24LST

    data Prcp /     0.0_RP, & ! 00LST
                    0.0_RP, & ! 01LST
                    0.0_RP, & ! 02LST
                    0.0_RP, & ! 03LST
                    0.1_RP, & ! 04LST
                    0.2_RP, & ! 05LST
                    3.0_RP, & ! 06LST
                    0.4_RP, & ! 07LST
                    0.0_RP, & ! 08LST
                    0.0_RP, & ! 09LST
                    0.0_RP, & ! 10LST
                    0.0_RP, & ! 11LST
                    0.0_RP, & ! 12LST
                    0.0_RP, & ! 13LST
                    0.0_RP, & ! 14LST
                    0.0_RP, & ! 15LST
                    0.0_RP, & ! 16LST
                    0.0_RP, & ! 17LST
                    0.0_RP, & ! 18LST
                    0.0_RP, & ! 19LST
                    0.0_RP, & ! 20LST
                    0.0_RP, & ! 21LST
                    0.0_RP, & ! 22LST
                    0.0_RP, & ! 23LST
                    0.0_RP  / ! 24LST

    real(RP) :: THETA
    real(RP) :: RovCP
    integer  :: i, j
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
                      + (        dsec ) * PT(tloc+1) )
          TMPA(i,j) = THETA * ( PRSA(i,j) / PRE00 )**RovCP  ! air temp, but now PRSA = 100000Pa

          UA  (i,j) = ( ( 1.0_RP-dsec ) * Wind(tloc  ) &
                      + (        dsec ) * Wind(tloc+1) )

          LWD (i,j) = RWD(i,j,I_LW,1) + RWD(i,j,I_LW,2)
          SWD (i,j) = RWD(i,j,I_SW,1) + RWD(i,j,I_SW,2)

          QVA (i,j) = ( ( 1.0_RP-dsec ) * Qvapor(tloc  ) &
                      + (        dsec ) * Qvapor(tloc+1) )
          QVA (i,j) = QVA(i,j) / (1.0_RP + QVA(i,j))         ! [mixing ratio->specific humidity]

          RAIN(i,j) = ( ( 1.0_RP-dsec ) * Prcp(tloc  ) &
                      + (        dsec ) * Prcp(tloc+1) ) / 3600.0_RP ! [mm/h->kg/m2/s]

!          if ( with_rain ) then
!             if (( tloc >= 1 ).and.( tloc < 10 )) then
!                RAIN(i,j) = 5.0_RP / 3600.0_RP
!             endif
!                RAIN(i,j) = 1.0_RP / 3600.0_RP
!          endif

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
