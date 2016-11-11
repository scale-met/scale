!-------------------------------------------------------------------------------
!> Module ideal topo
!!
!! @par Description
!!         This module prepares ideal topography for Dynamical Core Test
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_ideal_topo
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: IDEAL_topo

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: IDEAL_topo_JW
  private :: IDEAL_topo_Schar_Moderate
  private :: IDEAL_topo_Schar_Steep

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine IDEAL_topo( &
       lat, &
       lon, &
       Zsfc )
    use scale_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_KNONE, &
       ADM_lall,  &
       ADM_gall
    implicit none

    real(RP), intent(in)  :: lat (ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(in)  :: lon (ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(out) :: Zsfc(ADM_gall,ADM_KNONE,ADM_lall)

    character(len=H_SHORT) :: topo_type = ''

    namelist / IDEALTOPOPARAM / &
       topo_type

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ideal topo]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=IDEALTOPOPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** IDEALTOPOPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist IDEALTOPOPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=IDEALTOPOPARAM)

    if    ( topo_type == 'Schar_Moderate' ) then

       call IDEAL_topo_Schar_Moderate( lat (:,:,:), & ! [IN]
                                       lon (:,:,:), & ! [IN]
                                       Zsfc(:,:,:)  ) ! [OUT]

    elseif( topo_type == 'Schar_Steep' ) then

       call IDEAL_topo_Schar_Steep( lat (:,:,:), & ! [IN]
                                    lon (:,:,:), & ! [IN]
                                    Zsfc(:,:,:)  ) ! [OUT]

    elseif( topo_type == 'JW' ) then

       call IDEAL_topo_JW( lat (:,:,:), & ! [IN]
                           Zsfc(:,:,:)  ) ! [OUT]

    else
       write(*,*) 'xxx [IDEAL_topo] Not appropriate topo_type. STOP.'
       call PRC_MPIstop
    endif

    return
  end subroutine IDEAL_topo

  !-----------------------------------------------------------------------------
  !> Moderately-steep Schar-like circular mountain (Ref.: DCMIP 2012 eq.(48))
  subroutine IDEAL_topo_Schar_Moderate( &
       lat, &
       lon, &
       Zsfc )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       PI  => CONST_PI, &
       D2R => CONST_D2R
    use mod_adm, only: &
       ADM_KNONE, &
       ADM_lall,  &
       ADM_gall
    implicit none

    real(RP), intent(in)  :: lat (ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(in)  :: lon (ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(out) :: Zsfc(ADM_gall,ADM_KNONE,ADM_lall)

    real(RP) :: center_lon =  270.0_RP ! Longitude of Schar-type mountain center point [deg]
    real(RP) :: center_lat =    0.0_RP ! Latitude  of Schar-type mountain center point [deg]
    real(RP) :: H0         = 2000.0_RP ! Maximum Schar-type mountain height [m]
    real(RP) :: Rm_deg     =  135.0_RP ! Schar-type mountain radius     [deg]
    real(RP) :: QSIm_deg   =  11.25_RP ! Schar-type mountain wavelength [deg]

    namelist / IDEALTOPOPARAM_Schar_Moderate / &
       center_lon, &
       center_lat, &
       H0,         &
       Rm_deg,     &
       QSIm_deg

    real(RP) :: LAMBDA,  PHI
    real(RP) :: LAMBDAm, PHIm, Rm, QSIm
    real(RP) :: sinPHIm, cosPHIm
    real(RP) :: distance, mask

    integer :: g, l, K0
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[topo Schar Moderate]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=IDEALTOPOPARAM_Schar_Moderate,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** IDEALTOPOPARAM_Schar_Moderate is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist IDEALTOPOPARAM_Schar_Moderate. STOP.'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=IDEALTOPOPARAM_Schar_Moderate)

    K0 = ADM_KNONE

    LAMBDAm = center_lon * D2R ! [deg]->[rad]
    PHIm    = center_lat * D2R ! [deg]->[rad]
    Rm      = Rm_deg     * D2R ! [deg]->[rad]
    QSIm    = QSIm_deg   * D2R ! [deg]->[rad]
    sinPHIm = sin(PHIm)
    cosPHIm = cos(PHIm)

    do l = 1, ADM_lall
    do g = 1, ADM_gall
       LAMBDA = lon(g,K0,l)
       PHI    = lat(g,K0,l)

       distance = acos( sinPHIm * sin(PHI)                       &
                      + cosPHIm * cos(PHI) * cos(LAMBDA-LAMBDAm) )

       mask = 0.5_RP - sign(0.5_RP,distance-Rm) ! if distance > Rm, mask = 0

       Zsfc(g,ADM_KNONE,l) = H0/2.0_RP                              &
                           * ( 1.0_RP + cos( PI * distance / Rm ) ) &
                           * cos( PI * distance / QSIm )**2         &
                           * mask
    enddo
    enddo

    return
  end subroutine IDEAL_topo_Schar_Moderate

  !-----------------------------------------------------------------------------
  !> Steep Schar-type mountain (Ref.: DCMIP 2012 eq.(76))
  subroutine IDEAL_topo_Schar_Steep( &
       lat, &
       lon, &
       Zsfc )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       PI     => CONST_PI,     &
       D2R    => CONST_D2R,    &
       RADIUS => CONST_RADIUS
    use mod_adm, only: &
       ADM_KNONE, &
       ADM_lall,  &
       ADM_gall
    implicit none

    real(RP), intent(in)  :: lat (ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(in)  :: lon (ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(out) :: Zsfc(ADM_gall,ADM_KNONE,ADM_lall)

    real(RP) :: center_lon =   45.0_RP ! Longitude of Schar-type mountain center point [deg]
    real(RP) :: center_lat =    0.0_RP ! Latitude  of Schar-type mountain center point [deg]
    real(RP) :: H0         =  250.0_RP ! Maximum Schar-type mountain height [m]
    real(RP) :: d          = 5000.0_RP ! Schar-type mountain half-width [m]
    real(RP) :: QSI        = 4000.0_RP ! Schar-type mountain wavelength [m]

    namelist / IDEALTOPOPARAM_Schar_Steep / &
       center_lon, &
       center_lat, &
       H0,         &
       d,          &
       QSI

    real(RP) :: LAMBDA,  PHI
    real(RP) :: LAMBDAc, PHIc
    real(RP) :: sinPHIc, cosPHIc
    real(RP) :: distance

    integer :: g, l, K0
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[topo Schar Steep]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=IDEALTOPOPARAM_Schar_Steep,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** IDEALTOPOPARAM_Schar_Steep is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist IDEALTOPOPARAM_Schar_Steep. STOP.'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=IDEALTOPOPARAM_Schar_Steep)

    K0 = ADM_KNONE

    LAMBDAc = center_lon * D2R ! [deg]->[rad]
    PHIc    = center_lat * D2R ! [deg]->[rad]
    sinPHIc = sin(PHIc)
    cosPHIc = cos(PHIc)

    do l = 1, ADM_lall
    do g = 1, ADM_gall
       LAMBDA = lon(g,K0,l)
       PHI    = lat(g,K0,l)

       distance = RADIUS * acos( sinPHIc * sin(PHI)                       &
                               + cosPHIc * cos(PHI) * cos(LAMBDA-LAMBDAc) )

       Zsfc(g,ADM_KNONE,l) = H0                                  &
                           * exp( -(distance*distance) / (d*d) ) &
                           * cos( PI * distance / QSI )**2
    enddo
    enddo

    return
  end subroutine IDEAL_topo_Schar_Steep

  !-----------------------------------------------------------------------------
  !> mountain for JW06 testcase
  subroutine IDEAL_topo_JW( &
       lat, &
       Zsfc )
    use scale_const, only: &
       PI     => CONST_PI,     &
       GRAV   => CONST_GRAV,   &
       RADIUS => CONST_RADIUS, &
       OHM    => CONST_OHM
    use mod_adm, only: &
       ADM_KNONE, &
       ADM_lall,  &
       ADM_gall
    implicit none

    real(RP), intent(in)  :: lat (ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(out) :: Zsfc(ADM_gall,ADM_KNONE,ADM_lall)

    real(RP), parameter :: ETA0 = 0.252_RP ! Value of eta at a reference level (position of the jet)
    real(RP), parameter :: ETAs =   1.0_RP ! Value of eta at the surface
    real(RP), parameter :: u0   =  35.0_RP ! Maximum amplitude of the zonal wind

    real(RP) :: PHI, ETAv, u0cos32ETAv
    real(RP) :: f1, f2

    integer :: g, l, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    ETAv        = ( ETAs - ETA0 ) * PI/2.0_RP
    u0cos32ETAv = u0 * cos(ETAv)**(3.0_RP/2.0_RP)

    ! for globe
    do l = 1, ADM_lall
    do g = 1, ADM_gall
       PHI = lat(g,K0,l)

       f1 = -2.0_RP * sin(PHI)**6 * ( cos(PHI)**2 + 1.0_RP/3.0_RP ) + 10.0_RP/63.0_RP
       f2 = 8.0_RP/5.0_RP * cos(PHI)**3 * ( sin(PHI)**2 + 2.0_RP/3.0_RP ) - PI/4.0_RP

       Zsfc(g,k0,l) = u0cos32ETAv * ( u0cos32ETAv*f1 + RADIUS*OHM*f2 ) / GRAV
    enddo
    enddo

    return
  end subroutine IDEAL_topo_JW

end module mod_ideal_topo
