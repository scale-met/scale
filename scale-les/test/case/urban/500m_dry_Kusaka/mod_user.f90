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
  !logical, private :: with_rain = .false. !< add rain?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
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
    !if( IO_L ) write(IO_FID_LOG,*) '*** With rain? : ', with_rain

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
   use scale_const, only: &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, & 
       LH0   => CONST_LH0,   &    ! ELL : latent heat of vaporization [J/kg]
       Rvap  => CONST_Rvap        !< gas constant (water vapor) [J/kg/K]
    use scale_grid_real, only: &
       REAL_lon
    use mod_cpl_vars, only: &
       TMPA  => CPL_fromAtm_ATM_TEMP,    &
       QVA   => CPL_fromAtm_ATM_QV,      &
       UA    => CPL_fromAtm_ATM_U,       &
       VA    => CPL_fromAtm_ATM_V,       &
       WA    => CPL_fromAtm_ATM_W,       &
       RHOA  => CPL_fromAtm_ATM_DENS,    &
       PREC  => CPL_fromAtm_FLX_precip,  &
       SWD   => CPL_fromAtm_FLX_SW_dn,   &
       LWD   => CPL_fromAtm_FLX_LW_dn,   &
       PRES  => CPL_fromAtm_SFC_PRES,   &
       SHFLX => CPL_AtmUrb_ATM_FLX_SH,   &
       LHFLX => CPL_AtmUrb_ATM_FLX_LH,   &
       GHFLX => CPL_AtmUrb_URB_FLX_heat, &
       CNT_AtmUrb,                       &
       CNT_Urb
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

    real(RP) :: PT  (0:24)
    real(RP) :: Wind(0:24)
    real(RP) :: SW  (0:24)

    data SW / 0.0,0.0,0.0,0.0,0.0,0.0,50.0,240.0,420.0,600.0,690.0,765.0,800.0, &
              765.0,690.0,600.0,420.0,240.0,50.0,0.0,0.0,0.0,0.0,0.0,0.0/

    data PT / 28.0,27.8,27.65,27.5,27.35,27.2,27.1,27.5,27.85,28.25,28.8,29.4,30.0, &
              30.2,30.4,30.6,30.35,30.1,29.85,29.55,29.15,28.75,28.5,28.25,28.0/

    data Wind / 2.75,2.75,2.75,2.75,2.75,2.75,2.8,3.0,3.25,3.5,3.65,3.65,3.5, &
                3.4,3.27,3.15,3.05,2.95,2.85,2.8,2.75,2.7,2.72,2.75,2.75/

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( USER_do .AND. TIME_DOURBAN_step ) then

       CNT_AtmUrb = 0.0_RP
       CNT_Urb    = 0.0_RP

       QVA (:,:)  = 0.015_RP
       VA  (:,:)  = 0.0_RP
       WA  (:,:)  = 0.0_RP
       RHOA(:,:)  = 1.13_RP
       LWD (:,:)  = 400.0_RP
       PREC(:,:)  = 0.0_RP
       PRES(:,:)  = 101000.0_RP

       do j = 1, JA
       do i = 1, IA

          LON = REAL_lon(i,j) / D2R

          tloc = mod( int(NOWSEC/3600.0_RP)+int(LON/15.0_RP), 24 )

          dsec = mod(NOWSEC,3600.0_RP) / 3600.0_RP

          SWD (i,j) = ( ( 1.0_RP-dsec ) * SW(tloc  ) &
                      + (        dsec ) * SW(tloc+1) )
              
          TMPA(i,j) = ( ( 1.0_RP-dsec ) * PT(tloc  ) &
                      + (        dsec ) * PT(tloc+1) ) + TEM00

          UA  (i,j) = ( ( 1.0_RP-dsec ) * Wind(tloc  ) &
                      + (        dsec ) * Wind(tloc+1) )

       enddo
       enddo


       call HIST_in( TMPA (:,:), 'PT_urb',    'Potential temp',               'K',     dt_URB )
       call HIST_in( QVA  (:,:), 'QA_urb',    'Specific humidity',            'kg/kg', dt_URB )
       call HIST_in( UA   (:,:), 'UA_urb',    'Wind speed',                   'm/s',   dt_URB )
       call HIST_in( SWD  (:,:), 'SWD_urb',   'Downward shortwave radiation', 'W/m2',  dt_URB )
       call HIST_in( LWD  (:,:), 'LWD_urb',   'Downward longwave radiation',  'W/m2',  dt_URB )

       WORK(:,:) = PREC(:,:) * dt_URB
       call HIST_in( WORK (:,:), 'RAIN_urb',  'Precipitation',                'kg/m2', dt_URB )

       call HIST_in( SHFLX(:,:), 'SHFLX_urb', 'Sensible heat flux',           'W/m2',  dt_URB )
       call HIST_in( LHFLX(:,:), 'LHFLX_urb', 'Latent heat flux',             'W/m2',  dt_URB )

       WORK(:,:) = -GHFLX(:,:)
       call HIST_in( WORK(:,:),  'GHFLX_urb', 'Ground heat flux (down.pos.)', 'W/m2',  dt_URB )
    endif

    return
  end subroutine USER_step

end module mod_user
