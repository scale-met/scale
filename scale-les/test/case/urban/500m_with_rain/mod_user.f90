!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!      Put atmospheric data for urban test
!!      Test is based on Kusaka et al. (2000,BLM) 
!!
!! @author Team SCALE
!!
!!
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
  logical,  private, save :: USER_do = .false. !< do user step?

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
   use scale_const, only: &
       PI    => CONST_PI,     &
       TEM00 => CONST_TEM00
    use scale_grid_real, only: &
       XLON => REAL_lon, &
       XLAT => REAL_lat
    use mod_cpl_vars, only: &
       TMPA => CPL_fromAtm_ATM_TEMP,  &
       QVA  => CPL_fromAtm_ATM_QV,    &
       UA   => CPL_fromAtm_ATM_U,     &
       VA   => CPL_fromAtm_ATM_V,     &
       WA   => CPL_fromAtm_ATM_W,     &
       RHOA => CPL_fromAtm_ATM_DENS,   &
       PREC => CPL_fromAtm_FLX_precip, &
       SWD  => CPL_fromAtm_FLX_SW_dn,   &
       LWD  => CPL_fromAtm_FLX_LW_dn,   &
       SHFLX => CPL_AtmUrb_ATM_FLX_SH,  &
       LHFLX => CPL_AtmUrb_ATM_FLX_LH,  &
       GHFLX => CPL_AtmUrb_URB_FLX_heat
    use scale_time, only:   &
       TIME => TIME_NOWSEC, &   !< absolute sec
       TIME_DTSEC_URBAN !< time interval of urban step [sec]
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP) :: LON, LAT
    real(RP) :: tloc,dsec,DELT

    real(RP) :: PT(0:24)
    real(RP) :: Wind(0:24)
    real(RP) :: SW(0:24)
    data PT /28.0,27.8,27.65,27.5,27.35,27.2,27.1,27.5,27.85,28.25,28.8,29.4,30.0, &
             30.2,30.4,30.6,30.35,30.1,29.85,29.55,29.15,28.75,28.5,28.25,28.0/
    data Wind /2.75,2.75,2.75,2.75,2.75,2.75,2.8,3.0,3.25,3.5,3.65,3.65,3.5,       &
             3.4,3.27,3.15,3.05,2.95,2.85,2.8,2.75,2.7,2.72,2.75,2.75/
    data SW /0.0,0.0,0.0,0.0,0.0,0.0,50.0,240.0,420.0,600.0,690.0,765.0,800.0, &            
            765.0,690.0,600.0,420.0,240.0,50.0,0.0,0.0,0.0,0.0,0.0,0.0/

    integer :: i, j, k
    !---------------------------------------------------------------------------

    if ( USER_do ) then

       DELT=TIME_DTSEC_URBAN

       do j = 1, JA
       do i = 1, IA

         LON=XLON(i,j)/PI*180.0_RP
         tloc = mod( (int(TIME/(60.0_RP*60.0_RP)) + int(LON/15.0_RP)),24 )
         dsec = mod(TIME,3600.0_RP)

         TMPA(i,j) = (PT(tloc)*(3600.0_RP-dsec)+PT(tloc+1)*dsec)/3600.0_RP + TEM00
         QVA(i,j)  = 0.01_RP
         UA(i,j)   = (Wind(tloc)*(3600.0_RP-dsec)+Wind(tloc+1)*dsec)/3600.0_RP 
         VA(i,j)   = 0.0_RP
         WA(i,j)   = 0.0_RP
         RHOA(i,j) = 1.13_RP
 
         if (tloc < 10)then
            PREC(i,j) = 5.0_RP/3600.0_RP
         else
            PREC(i,j) = 0.0_RP/3600.0_RP
         endif

         SWD(i,j)  = (SW(tloc)*(3600.0_RP-dsec)+SW(tloc+1)*dsec)/3600.0_RP 
         LWD(i,j)  = 400.0_RP

       !print *,'user',TMPA(i,j),QVA(i,j),SWD(i,j)

       enddo
       enddo

       call HIST_in( SWD(:,:),  'SWD_urb', 'Downward shortwave radiation', 'W/m2', TIME_DTSEC_URBAN)
       call HIST_in( LWD(:,:),  'LWD_urb', 'Downward longwave radiation',  'W/m2', TIME_DTSEC_URBAN)
       call HIST_in( TMPA(:,:), 'PT_urb',  'Potential temp', 'K',   TIME_DTSEC_URBAN)
       call HIST_in( UA(:,:),   'UA_urb',  'Wind speed',     'm/s', TIME_DTSEC_URBAN)

       call HIST_in( SHFLX(:,:), 'SHFLX_urb', 'Sensible heat flux', 'W/m2', TIME_DTSEC_URBAN)
       call HIST_in( LHFLX(:,:), 'LHFLX_urb', 'Latent heat flux',   'W/m2', TIME_DTSEC_URBAN)
       call HIST_in( GHFLX(:,:), 'GHFLX_urb', 'Ground heat flux',   'W/m2', TIME_DTSEC_URBAN)

    endif

    return
  end subroutine USER_step

end module mod_user
