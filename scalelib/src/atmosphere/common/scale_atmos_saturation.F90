!-------------------------------------------------------------------------------
!> module atmosphere / saturation
!!
!! @par Description
!!          Saturation module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_saturation
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_const, only: &
     Rvap   => CONST_Rvap,   &
     PSAT0  => CONST_PSAT0,  &
     EPSvap => CONST_EPSvap, &
     TEM00  => CONST_TEM00
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_SATURATION_setup

  public :: ATMOS_SATURATION_alpha

  public :: ATMOS_SATURATION_psat_all
  public :: ATMOS_SATURATION_psat_liq
  public :: ATMOS_SATURATION_psat_ice

  public :: ATMOS_SATURATION_psat2qsat_pres
  public :: ATMOS_SATURATION_psat2qsat_dens

  public :: ATMOS_SATURATION_pres2qsat_all
  public :: ATMOS_SATURATION_pres2qsat_liq
  public :: ATMOS_SATURATION_pres2qsat_ice

  public :: ATMOS_SATURATION_dens2qsat_all
  public :: ATMOS_SATURATION_dens2qsat_liq
  public :: ATMOS_SATURATION_dens2qsat_ice

  public :: ATMOS_SATURATION_dalphadT

  public :: ATMOS_SATURATION_dqs_dtem_dens_liq
  public :: ATMOS_SATURATION_dqs_dtem_dens_ice
  public :: ATMOS_SATURATION_dqs_dtem_dpre_liq
  public :: ATMOS_SATURATION_dqs_dtem_dpre_ice

  public :: ATMOS_SATURATION_tdew_liq

  public :: ATMOS_SATURATION_pote

  public :: ATMOS_SATURATION_moist_conversion_dens_liq
  public :: ATMOS_SATURATION_moist_conversion_dens_all
  public :: ATMOS_SATURATION_moist_conversion_pres_liq

  interface ATMOS_SATURATION_alpha
     module procedure ATMOS_SATURATION_alpha_0D
     module procedure ATMOS_SATURATION_alpha_1D
     module procedure ATMOS_SATURATION_alpha_3D
  end interface ATMOS_SATURATION_alpha

  interface ATMOS_SATURATION_psat_all
     module procedure ATMOS_SATURATION_psat_all_0D
     module procedure ATMOS_SATURATION_psat_all_1D
     module procedure ATMOS_SATURATION_psat_all_2D
     module procedure ATMOS_SATURATION_psat_all_3D
  end interface ATMOS_SATURATION_psat_all
  interface ATMOS_SATURATION_psat_liq
     module procedure ATMOS_SATURATION_psat_liq_0D
     module procedure ATMOS_SATURATION_psat_liq_1D
     module procedure ATMOS_SATURATION_psat_liq_2D
     module procedure ATMOS_SATURATION_psat_liq_3D
  end interface ATMOS_SATURATION_psat_liq
  interface ATMOS_SATURATION_psat_ice
     module procedure ATMOS_SATURATION_psat_ice_0D
     module procedure ATMOS_SATURATION_psat_ice_1D
     module procedure ATMOS_SATURATION_psat_ice_2D
     module procedure ATMOS_SATURATION_psat_ice_3D
  end interface ATMOS_SATURATION_psat_ice

  interface ATMOS_SATURATION_psat2qsat_pres
     module procedure ATMOS_SATURATION_psat2qsat_pres_0D
  end interface ATMOS_SATURATION_psat2qsat_pres
  interface ATMOS_SATURATION_psat2qsat_dens
     module procedure ATMOS_SATURATION_psat2qsat_dens_0D
  end interface ATMOS_SATURATION_psat2qsat_dens

  interface ATMOS_SATURATION_pres2qsat_all
     module procedure ATMOS_SATURATION_pres2qsat_all_0D
     module procedure ATMOS_SATURATION_pres2qsat_all_1D
     module procedure ATMOS_SATURATION_pres2qsat_all_2D
     module procedure ATMOS_SATURATION_pres2qsat_all_3D
  end interface ATMOS_SATURATION_pres2qsat_all
  interface ATMOS_SATURATION_pres2qsat_liq
     module procedure ATMOS_SATURATION_pres2qsat_liq_0D
     module procedure ATMOS_SATURATION_pres2qsat_liq_1D
     module procedure ATMOS_SATURATION_pres2qsat_liq_3D
  end interface ATMOS_SATURATION_pres2qsat_liq
  interface ATMOS_SATURATION_pres2qsat_ice
     module procedure ATMOS_SATURATION_pres2qsat_ice_0D
     module procedure ATMOS_SATURATION_pres2qsat_ice_1D
     module procedure ATMOS_SATURATION_pres2qsat_ice_3D
  end interface ATMOS_SATURATION_pres2qsat_ice

  interface ATMOS_SATURATION_dens2qsat_all
     module procedure ATMOS_SATURATION_dens2qsat_all_0D
     module procedure ATMOS_SATURATION_dens2qsat_all_1D
     module procedure ATMOS_SATURATION_dens2qsat_all_3D
  end interface ATMOS_SATURATION_dens2qsat_all
  interface ATMOS_SATURATION_dens2qsat_liq
     module procedure ATMOS_SATURATION_dens2qsat_liq_0D
     module procedure ATMOS_SATURATION_dens2qsat_liq_1D
     module procedure ATMOS_SATURATION_dens2qsat_liq_3D
  end interface ATMOS_SATURATION_dens2qsat_liq
  interface ATMOS_SATURATION_dens2qsat_ice
     module procedure ATMOS_SATURATION_dens2qsat_ice_0D
     module procedure ATMOS_SATURATION_dens2qsat_ice_1D
     module procedure ATMOS_SATURATION_dens2qsat_ice_3D
  end interface ATMOS_SATURATION_dens2qsat_ice

  interface ATMOS_SATURATION_dalphadT
     module procedure ATMOS_SATURATION_dalphadT_0D
     module procedure ATMOS_SATURATION_dalphadT_1D
     module procedure ATMOS_SATURATION_dalphadT_3D
  end interface ATMOS_SATURATION_dalphadT

  interface ATMOS_SATURATION_dqs_dtem_dens_liq
     module procedure ATMOS_SATURATION_dqs_dtem_dens_liq_0D
     module procedure ATMOS_SATURATION_dqs_dtem_dens_liq_1D
     module procedure ATMOS_SATURATION_dqs_dtem_dens_liq_3D
  end interface ATMOS_SATURATION_dqs_dtem_dens_liq
  interface ATMOS_SATURATION_dqs_dtem_dens_ice
     module procedure ATMOS_SATURATION_dqs_dtem_dens_ice_0D
     module procedure ATMOS_SATURATION_dqs_dtem_dens_ice_1D
     module procedure ATMOS_SATURATION_dqs_dtem_dens_ice_3D
  end interface ATMOS_SATURATION_dqs_dtem_dens_ice
  interface ATMOS_SATURATION_dqs_dtem_dpre_liq
     module procedure ATMOS_SATURATION_dqs_dtem_dpre_liq_0D
     module procedure ATMOS_SATURATION_dqs_dtem_dpre_liq_1D
     module procedure ATMOS_SATURATION_dqs_dtem_dpre_liq_3D
  end interface ATMOS_SATURATION_dqs_dtem_dpre_liq
  interface ATMOS_SATURATION_dqs_dtem_dpre_ice
     module procedure ATMOS_SATURATION_dqs_dtem_dpre_ice_0D
     module procedure ATMOS_SATURATION_dqs_dtem_dpre_ice_1D
     module procedure ATMOS_SATURATION_dqs_dtem_dpre_ice_3D
  end interface ATMOS_SATURATION_dqs_dtem_dpre_ice

  interface ATMOS_SATURATION_tdew_liq
     module procedure ATMOS_SATURATION_tdew_liq_0D
     module procedure ATMOS_SATURATION_tdew_liq_1D
     module procedure ATMOS_SATURATION_tdew_liq_3D
  end interface ATMOS_SATURATION_tdew_liq

  interface ATMOS_SATURATION_pote
     module procedure ATMOS_SATURATION_pote_0D
     module procedure ATMOS_SATURATION_pote_1D
     module procedure ATMOS_SATURATION_pote_3D
  end interface ATMOS_SATURATION_pote

  interface ATMOS_SATURATION_moist_conversion_dens_liq
     module procedure ATMOS_SATURATION_moist_conversion_dens_liq_0D
  end interface ATMOS_SATURATION_moist_conversion_dens_liq
  interface ATMOS_SATURATION_moist_conversion_dens_all
     module procedure ATMOS_SATURATION_moist_conversion_dens_all_0D
  end interface ATMOS_SATURATION_moist_conversion_dens_all
  interface ATMOS_SATURATION_moist_conversion_pres_liq
     module procedure ATMOS_SATURATION_moist_conversion_pres_liq_0D
  end interface ATMOS_SATURATION_moist_conversion_pres_liq

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
  real(RP), private, parameter :: TEM_MIN = 10.0_RP !> calculation of dew point is ommited under this temperature

  real(RP), private :: ATMOS_SATURATION_ULIMIT_TEMP = 273.15_RP !< upper limit temperature
  real(RP), private :: ATMOS_SATURATION_LLIMIT_TEMP = 233.15_RP !< lower limit temperature

  real(RP), private :: RTEM00         !> inverse of TEM00
  real(RP), private :: dalphadT_const !> d(alfa)/dt
  real(RP), private :: psat_min_liq   !> psat_liq for TEM_MIN
  real(RP), private :: psat_min_ice   !> psat_ice for TEM_MIN

  real(RP), private :: CPovR_liq
  real(RP), private :: CPovR_ice
  real(RP), private :: CVovR_liq
  real(RP), private :: CVovR_ice
  real(RP), private :: LovR_liq
  real(RP), private :: LovR_ice

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_SATURATION_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       CPvap => CONST_CPvap, &
       CVvap => CONST_CVvap, &
       CL    => CONST_CL,    &
       CI    => CONST_CI,    &
       LHV00 => CONST_LHV00, &
       LHS00 => CONST_LHS00, &
       LHV0  => CONST_LHV0,  &
       LHS0  => CONST_LHS0,  &
       CONST_THERMODYN_TYPE
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_setup
    implicit none

    namelist / PARAM_ATMOS_SATURATION / &
       ATMOS_SATURATION_ULIMIT_TEMP, &
       ATMOS_SATURATION_LLIMIT_TEMP

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_SATURATION_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_SATURATION,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_SATURATION_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_SATURATION_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_SATURATION. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_SATURATION)

    RTEM00 = 1.0_RP / TEM00

    if ( CONST_THERMODYN_TYPE == 'EXACT' ) then

       CPovR_liq = ( CPvap - CL ) / Rvap
       CPovR_ice = ( CPvap - CI ) / Rvap
       CVovR_liq = ( CVvap - CL ) / Rvap
       CVovR_ice = ( CVvap - CI ) / Rvap

       LovR_liq  = LHV00 / Rvap
       LovR_ice  = LHS00 / Rvap

    elseif( CONST_THERMODYN_TYPE == 'SIMPLE' ) then

       CPovR_liq = 0.0_RP
       CPovR_ice = 0.0_RP
       CVovR_liq = 0.0_RP
       CVovR_ice = 0.0_RP

       LovR_liq  = LHV0 / Rvap
       LovR_ice  = LHS0 / Rvap

    endif

    dalphadT_const = 1.0_RP / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )

    LOG_NEWLINE
    LOG_INFO("ATMOS_SATURATION_setup",'(1x,A,F7.2,A,F7.2)') 'Temperature range for liquid/ice mixture : ', &
                                                      ATMOS_SATURATION_LLIMIT_TEMP, ' - ', &
                                                      ATMOS_SATURATION_ULIMIT_TEMP

    call ATMOS_HYDROMETEOR_setup

    call ATMOS_SATURATION_psat_liq( TEM_MIN, psat_min_liq ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat_ice( TEM_MIN, psat_min_ice ) ! [IN], [OUT]

    return
  end subroutine ATMOS_SATURATION_setup

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (0D)
  subroutine ATMOS_SATURATION_alpha_0D( &
       temp, &
       alpha )
    implicit none

    real(RP), intent(in)  :: temp  !< temperature [K]
    real(RP), intent(out) :: alpha !< liquid/ice separation factor (0-1)
    !---------------------------------------------------------------------------

    alpha = ( temp                         - ATMOS_SATURATION_LLIMIT_TEMP ) &
          / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )

    alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

    return
  end subroutine ATMOS_SATURATION_alpha_0D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_alpha_1D( &
       KA, KS, KE, &
       temp, &
       alpha )
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: temp (KA) !< temperature [K]

    real(RP), intent(out) :: alpha(KA) !< liquid/ice separation factor (0-1)

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_alpha_0D( temp(k), & ! [IN]
                                       alpha(k) ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_alpha_1D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (3D)
  subroutine ATMOS_SATURATION_alpha_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, &
       alpha )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp (KA,IA,JA) !< temperature [K]
    real(RP), intent(out) :: alpha(KA,IA,JA) !< liquid/ice separation factor (0-1)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_alpha_0D( temp(k,i,j), & ! [IN]
                                       alpha(k,i,j) ) ! [OUT]

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_alpha_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (0D)
  subroutine ATMOS_SATURATION_psat_all_0D( &
       temp, &
       psat  )
    implicit none

    real(RP), intent(in)  :: temp !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]

    real(RP) :: alpha, psatl, psati
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_alpha   ( temp, alpha )
    call ATMOS_SATURATION_psat_liq( temp, psatl )
    call ATMOS_SATURATION_psat_ice( temp, psati )

    psat = psatl * (          alpha ) &
         + psati * ( 1.0_RP - alpha )

    return
  end subroutine ATMOS_SATURATION_psat_all_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_psat_all_1D( &
       KA, KS, KE, &
       temp, &
       psat  )
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA) !< temperature               [K]
    real(RP), intent(out) :: psat(KA) !< saturation vapor pressure [Pa]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_psat_all_0D( temp(k), & ! [IN]
                                          psat(k)  ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_all_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (2D)
  subroutine ATMOS_SATURATION_psat_all_2D( &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       psat  )
    implicit none
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(IA,JA) !< temperature               [K]

    real(RP), intent(out) :: psat(IA,JA) !< saturation vapor pressure [Pa]

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_SATURATION_psat_all_0D( temp(i,j), & ! [IN]
                                          psat(i,j)  ) ! [OUT]
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_all_2D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (3D)
  subroutine ATMOS_SATURATION_psat_all_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       psat  )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature               [K]
    real(RP), intent(out) :: psat(KA,IA,JA) !< saturation vapor pressure [Pa]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_psat_all_0D( temp(k,i,j), & ! [IN]
                                          psat(k,i,j)  ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_all_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (0D)
  subroutine ATMOS_SATURATION_psat_liq_0D( &
       temp, &
       psat  )
    implicit none

    real(RP), intent(in)  :: temp !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( temp * RTEM00 )**CPovR_liq             &
                 * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp ) )

    return
  end subroutine ATMOS_SATURATION_psat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_psat_liq_1D( &
       KA, KS, KE, &
       temp, &
       psat  )
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA) !< temperature               [K]
    real(RP), intent(out) :: psat(KA) !< saturation vapor pressure [Pa]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_psat_liq_0D( temp(k), &
                                          psat(k) )
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_liq_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (2D)
  subroutine ATMOS_SATURATION_psat_liq_2D( &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       psat  )
    implicit none

    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(IA,JA) !< temperature               [K]
    real(RP), intent(out) :: psat(IA,JA) !< saturation vapor pressure [Pa]

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_SATURATION_psat_liq_0D( temp(i,j), & ! [IN]
                                          psat(i,j)  ) ! [OUT]
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_liq_2D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (3D)
  subroutine ATMOS_SATURATION_psat_liq_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       psat  )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature               [K]
    real(RP), intent(out) :: psat(KA,IA,JA) !< saturation vapor pressure [Pa]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_psat_liq_0D( temp(k,i,j), & ! [IN]
                                          psat(k,i,j)  ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_liq_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (0D)
  subroutine ATMOS_SATURATION_psat_ice_0D( &
       temp, &
       psat  )
    implicit none

    real(RP), intent(in)  :: temp !< temperature               [K]

    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( temp * RTEM00 )**CPovR_ice             &
                 * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp ) )

    return
  end subroutine ATMOS_SATURATION_psat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_psat_ice_1D( &
       KA, KS, KE, &
       temp, &
       psat  )
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA) !< temperature               [K]
    real(RP), intent(out) :: psat(KA) !< saturation vapor pressure [Pa]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_psat_ice_0D( temp(k), & ! [IN]
                                          psat(k)  ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_ice_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (2D)
  subroutine ATMOS_SATURATION_psat_ice_2D( &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       psat  )
    implicit none

    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(IA,JA) !< temperature               [K]
    real(RP), intent(out) :: psat(IA,JA) !< saturation vapor pressure [Pa]

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_SATURATION_psat_ice_0D( temp(i,j), & ! [IN]
                                          psat(i,j)  ) ! [OUT]
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_ice_2D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (3D)
  subroutine ATMOS_SATURATION_psat_ice_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       psat  )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature               [K]
    real(RP), intent(out) :: psat(KA,IA,JA) !< saturation vapor pressure [Pa]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_psat_ice_0D( temp(k,i,j), & ! [IN]
                                          psat(k,i,j)  ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_ice_3D


  !-----------------------------------------------------------------------------
  !> calc saturation pressure -> saturation vapor mass under constant pressure
  subroutine ATMOS_SATURATION_psat2qsat_pres_0D( &
       psat, pres, qdry, &
       qsat              )
    implicit none

    real(RP), intent(in)  :: psat !< saturation pressure [Pa]
    real(RP), intent(in)  :: pres !< pressure            [Pa]
    real(RP), intent(in)  :: qdry !< dry air mass ratio  [kg/kg]

    real(RP), intent(out) :: qsat !< saturation vapor mass [kg/kg]

    ! ! qdry is assumed to be 1 - qsat
    ! qsat = EPSvap * psat / ( pres - ( 1.0_RP-EPSvap ) * psat )

    qsat = EPSvap * qdry * psat / ( pres - psat )

    return
  end subroutine ATMOS_SATURATION_psat2qsat_pres_0D
  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,0D)
  subroutine ATMOS_SATURATION_pres2qsat_all_0D( &
       temp, pres, qdry, &
       qsat              )
    implicit none

    real(RP), intent(in)  :: temp !< temperature           [K]
    real(RP), intent(in)  :: pres !< pressure              [Pa]
    real(RP), intent(in)  :: qdry !< dry air mass ratio    [kg/kg]

    real(RP), intent(out) :: qsat !< saturation vapor mass [kg/kg]

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_all_0D( temp, psat ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat2qsat_pres_0D( psat, pres, qdry, & ! [IN]
                                             qsat              ) ! [OUT]

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_0D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_pres2qsat_all_1D( &
       KA, KS, KE, &
       temp, pres, qdry, &
       qsat              )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA) !< temperature           [K]
    real(RP), intent(in)  :: pres(KA) !< pressure              [Pa]
    real(RP), intent(in)  :: qdry(KA) !< dry air mass ratio    [kg/kg]

    real(RP), intent(out) :: qsat(KA) !< saturation vapor mass [kg/kg]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_pres2qsat_all_0D( temp(k), pres(k), qdry(k), & ! [IN]
                                               qsat(k)                    ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_1D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,2D)
  subroutine ATMOS_SATURATION_pres2qsat_all_2D( &
       IA, IS, IE, JA, JS, JE, &
       temp, pres, qdry, &
       qsat              )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp(IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: pres(IA,JA) !< pressure              [Pa]
    real(RP), intent(in)  :: qdry(IA,JA) !< dry air mass ratio    [kg/kg]

    real(RP), intent(out) :: qsat(IA,JA) !< saturation vapor mass [kg/kg]

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_SATURATION_pres2qsat_all_0D( temp(i,j), pres(i,j), qdry(i,j), &
                                               qsat(i,j)                        )
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_2D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,3D)
  subroutine ATMOS_SATURATION_pres2qsat_all_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, pres, qdry, &
       qsat        )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)  :: qdry(KA,IA,JA) !< dry air mass ratio    [kg/kg]

    real(RP), intent(out) :: qsat(KA,IA,JA) !< saturation vapor mass [kg/kg]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_pres2qsat_all_0D( temp(k,i,j), pres(k,i,j), qdry(k,i,j), & ! [IN]
                                               qsat(k,i,j)                            ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_3D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid,0D)
  subroutine ATMOS_SATURATION_pres2qsat_liq_0D( &
       temp, pres, qdry, &
       qsat              )
    implicit none

    real(RP), intent(in)  :: temp !< temperature           [K]
    real(RP), intent(in)  :: pres !< pressure              [Pa]
    real(RP), intent(in)  :: qdry !< dry air mass ratio    [kg/kg]

    real(RP), intent(out) :: qsat !< saturation vapor mass [kg/kg]

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_liq( temp, psat ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat2qsat_pres( psat, pres, qdry, & ! [IN]
                                          qsat              ) ! [OUT]

    return
  end subroutine ATMOS_SATURATION_pres2qsat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid,1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_pres2qsat_liq_1D( &
       KA, KS, KE, &
       temp, pres, qdry, &
       qsat              )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA) !< temperature           [K]
    real(RP), intent(in)  :: pres(KA) !< pressure              [Pa]
    real(RP), intent(in)  :: qdry(KA) !< dry air mass ratio    [kg/kg]

    real(RP), intent(out) :: qsat(KA) !< saturation vapor mass [kg/kg]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_pres2qsat_liq_0D( temp(k), pres(k), qdry(k), & ! [IN]
                                               qsat(k)                    ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_liq_1D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid,3D)
  subroutine ATMOS_SATURATION_pres2qsat_liq_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, pres, qdry, &
       qsat              )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)  :: qdry(KA,IA,JA) !< dry air mass ratio    [kg/kg]

    real(RP), intent(out) :: qsat(KA,IA,JA) !< saturation vapor mass [kg/kg]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_pres2qsat_liq_0D( temp(k,i,j), pres(k,i,j), qdry(k,i,j), & ! [IN]
                                               qsat(k,i,j)                            ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_liq_3D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (ice,0D)
  subroutine ATMOS_SATURATION_pres2qsat_ice_0D( &
       temp, pres, qdry, &
       qsat              )
    implicit none

    real(RP), intent(in)  :: temp !< temperature           [K]
    real(RP), intent(in)  :: pres !< pressure              [Pa]
    real(RP), intent(in)  :: qdry !< dry air mass ratio    [kg/kg]

    real(RP), intent(out) :: qsat !< saturation vapor mass [kg/kg]

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_ice( temp, psat ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat2qsat_pres( psat, pres, qdry, & ! [IN]
                                          qsat              ) ! [OUT]

    return
  end subroutine ATMOS_SATURATION_pres2qsat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (ice,1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_pres2qsat_ice_1D( &
       KA, KS, KE, &
       temp, pres, qdry, &
       qsat              )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(in)  :: pres(KA)
    real(RP), intent(in)  :: qdry(KA)

    real(RP), intent(out) :: qsat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_pres2qsat_ice_0D( temp(k), pres(k), qdry(k), & ! [IN]
                                               qsat(k)                    ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_ice_1D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (ice,3D)
  subroutine ATMOS_SATURATION_pres2qsat_ice_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, pres, qdry, &
       qsat        )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: pres(KA,IA,JA)
    real(RP), intent(in)  :: qdry(KA,IA,JA)

    real(RP), intent(out) :: qsat(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_pres2qsat_ice_0D( temp(k,i,j), pres(k,i,j), qdry(k,i,j), & ! [IN]
                                               qsat(k,i,j)                            ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_ice_3D

  !-----------------------------------------------------------------------------
  !> calc saturation pressure -> saturation vapor mass under constant density (volume)
  subroutine ATMOS_SATURATION_psat2qsat_dens_0D( &
       psat, temp, dens, &
       qsat              )
    implicit none
    real(RP), intent(in)  :: psat
    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP), intent(out) :: qsat

    qsat = psat / ( dens * Rvap * temp )

    return
  end subroutine ATMOS_SATURATION_psat2qsat_dens_0D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid/ice mixture,0D)
  subroutine ATMOS_SATURATION_dens2qsat_all_0D( &
       temp, &
       dens, &
       qsat  )
    implicit none
    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP), intent(out) :: qsat

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_all( temp, psat ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat2qsat_dens_0D( psat, temp, dens, & ! [IN]
                                             qsat              ) ! [OUT]

    return
  end subroutine ATMOS_SATURATION_dens2qsat_all_0D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid/ice mixture,1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_dens2qsat_all_1D( &
       KA, KS, KE, &
       temp, dens, &
       qsat        )
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(in)  :: dens(KA)

    real(RP), intent(out) :: qsat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_dens2qsat_all_0D( temp(k), dens(k), & ! [IN]
                                               qsat(k)           ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_all_1D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid/ice mixture,3D)
  subroutine ATMOS_SATURATION_dens2qsat_all_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       dens, &
       qsat  )
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: dens(KA,IA,JA)

    real(RP), intent(out) :: qsat(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_dens2qsat_all_0D( temp(k,i,j), dens(k,i,j), & ! [IN]
                                               qsat(k,i,j)               ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_all_3D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid,0D)
  subroutine ATMOS_SATURATION_dens2qsat_liq_0D( &
       temp, &
       dens, &
       qsat  )
    implicit none
    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP), intent(out) :: qsat

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_liq_0D( temp, psat ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat2qsat_dens_0D( psat, temp, dens, & ! [IN]
                                             qsat              ) ! [OUT]

    return
  end subroutine ATMOS_SATURATION_dens2qsat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid,1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_dens2qsat_liq_1D( &
       KA, KS, KE, &
       temp, dens, &
       qsat        )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(in)  :: dens(KA)

    real(RP), intent(out) :: qsat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_dens2qsat_liq_0D( temp(k), dens(k), & ! [IN]
                                               qsat(k)           ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_liq_1D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid,3D)
  subroutine ATMOS_SATURATION_dens2qsat_liq_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, dens, &
       qsat        )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: dens(KA,IA,JA)

    real(RP), intent(out) :: qsat(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_dens2qsat_liq_0D( temp(k,i,j), dens(k,i,j), & ! [IN]
                                               qsat(k,i,j)               ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_liq_3D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (ice,0D)
  subroutine ATMOS_SATURATION_dens2qsat_ice_0D( &
       temp, &
       dens, &
       qsat  )
    implicit none
    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP), intent(out) :: qsat

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_ice( temp, psat ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat2qsat_dens( psat, temp, dens, & ! [IN]
                                          qsat              ) ! [OUT]
    return
  end subroutine ATMOS_SATURATION_dens2qsat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (ice,1D)
!OCL SERIAL
  subroutine ATMOS_SATURATION_dens2qsat_ice_1D( &
       KA, KS, KE, &
       temp, dens, &
       qsat        )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(in)  :: dens(KA)

    real(RP), intent(out) :: qsat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_dens2qsat_ice_0D( temp(k), dens(k), & ! [IN]
                                               qsat(k)           ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_ice_1D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (ice,3D)
  subroutine ATMOS_SATURATION_dens2qsat_ice_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, dens, &
       qsat        )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: dens(KA,IA,JA)

    real(RP), intent(out) :: qsat(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_dens2qsat_ice_0D( temp(k,i,j), dens(k,i,j), & ! [IN]
                                               qsat(k,i,j)               ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_ice_3D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(temp), 0D
  subroutine ATMOS_SATURATION_dalphadT_0D( &
       temp,     &
       dalpha_dT )
    implicit none
    real(RP), intent(in)  :: temp

    real(RP), intent(out) :: dalpha_dT

    real(RP) :: lim1, lim2
    !---------------------------------------------------------------------------

    ! if Tup < temp, dalpha/dT = 0 (no slope)
    lim1 = 0.5_RP + sign( 0.5_RP, ATMOS_SATURATION_ULIMIT_TEMP - temp )
    ! if Tdn > temp, dalpha/dT = 0 (no slope)
    lim2 = 0.5_RP + sign( 0.5_RP, temp - ATMOS_SATURATION_LLIMIT_TEMP )

    dalpha_dT = dalphadT_const * lim1 * lim2

    return
  end subroutine ATMOS_SATURATION_dalphadT_0D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(temp), 1D
!OCL SERIAL
  subroutine ATMOS_SATURATION_dalphadT_1D( &
       KA, KS, KE, &
       temp,     &
       dalpha_dT )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp     (KA)

    real(RP), intent(out) :: dalpha_dT(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_dalphadT_0D( temp(k), dalpha_dT(k) ) ! [IN], [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_dalphadT_1D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(temp), 3D
  subroutine ATMOS_SATURATION_dalphadT_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp,     &
       dalpha_dT )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp     (KA,IA,JA)

    real(RP), intent(out) :: dalpha_dT(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_dalphadT_0D( temp(k,i,j), dalpha_dT(k,i,j) ) ! [IN], [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dalphadT_3D

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water at constant density
  subroutine ATMOS_SATURATION_dqs_dtem_dens_liq_0D( &
       temp, dens, &
       dqsdtem, qsat )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    implicit none

    real(RP), intent(in)  :: temp !< temperature           [K]
    real(RP), intent(in)  :: dens !< temperature           [K]

    real(RP), intent(out) :: dqsdtem !< (d qsw/d T)_{rho}
    real(RP), intent(out), optional :: qsat

    real(RP) :: LHV ! latent heat of vaporization [J/kg]
    real(RP) :: psat
    !---------------------------------------------------------------------------

    call HYDROMETEOR_LHV( temp, LHV )
    call ATMOS_SATURATION_psat_liq( temp, psat ) ! [IN], [OUT]

    dqsdtem = psat / ( dens* Rvap * temp**2 ) &
            * ( LHV / ( Rvap * temp ) - 1.0_RP )

    if ( present(qsat) ) &
         call ATMOS_SATURATION_psat2qsat_dens( psat, temp, dens, & ! [IN]
                                               qsat              ) ! [OUT]

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dens_liq_0D

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_SATURATION_dqs_dtem_dens_liq_1D( &
       KA, KS, KE, &
       temp, dens, &
       dqsdtem     )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp   (KA) !< temperature           [K]
    real(RP), intent(in)  :: dens   (KA) !< temperature           [K]

    real(RP), intent(out) :: dqsdtem(KA) !< (d qsw/d T)_{rho}

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_dqs_dtem_dens_liq_0D( temp(k), dens(k), & ! [IN]
                                                   dqsdtem(k)            ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dens_liq_1D

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water at constant density
  subroutine ATMOS_SATURATION_dqs_dtem_dens_liq_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, dens, &
       dqsdtem     )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp   (KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: dens   (KA,IA,JA) !< temperature           [K]

    real(RP), intent(out) :: dqsdtem(KA,IA,JA) !< (d qsw/d T)_{rho}

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_dqs_dtem_dens_liq_0D( temp(k,i,j), dens(k,i,j), & ! [IN]
                                                   dqsdtem(k,i,j)            ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dens_liq_3D

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice at constant density
  subroutine ATMOS_SATURATION_dqs_dtem_dens_ice_0D( &
       temp, dens,   &
       dqsdtem, qsat )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHS => ATMOS_HYDROMETEOR_LHS
    implicit none

    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP), intent(out) :: dqsdtem
    real(RP), intent(out), optional :: qsat

    real(RP) :: LHS ! latent heat of sublimation  [J/kg]
    real(RP) :: psat

    !---------------------------------------------------------------------------

    call HYDROMETEOR_LHS( temp, LHS ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat_ice( temp, psat ) ! [IN], [OUT]

    dqsdtem = psat / ( dens * Rvap * temp**2 ) &
            * ( LHS / ( Rvap * temp ) - 1.0_RP )

    if ( present(qsat) ) &
         call ATMOS_SATURATION_psat2qsat_dens( psat, temp, dens, & ! [IN]
                                               qsat              ) ! [OUT]

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dens_ice_0D

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_SATURATION_dqs_dtem_dens_ice_1D( &
       KA, KS, KE, &
       temp, dens, &
       dqsdtem     )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp   (KA)
    real(RP), intent(in)  :: dens   (KA)

    real(RP), intent(out) :: dqsdtem(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_dqs_dtem_dens_ice_0D( temp(k), dens(k), & ! [IN]
                                                   dqsdtem(k)        ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dens_ice_1D

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice at constant density
  subroutine ATMOS_SATURATION_dqs_dtem_dens_ice_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, dens, &
       dqsdtem     )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: dens   (KA,IA,JA)

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_dqs_dtem_dens_ice_0D( temp(k,i,j), dens(k,i,j), & ! [IN]
                                                   dqsdtem(k,i,j)            ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dens_ice_3D

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_all at constant density
  subroutine ATMOS_SATURATION_dqs_dtem_dens_all_0D( &
       temp, dens,               &
       dqsat_dT,                 &
       qsat, qsat_liq, qsat_ice, &
       alpha                     )
    implicit none

    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP), intent(out) :: dqsat_dT

    real(RP), intent(out), optional :: qsat
    real(RP), intent(out), optional :: qsat_liq
    real(RP), intent(out), optional :: qsat_ice
    real(RP), intent(out), optional :: alpha

    real(RP) :: qsat_liq_, qsat_ice_, alpha_
    real(RP) :: dqsat_dT_liq, dqsat_dT_ice, dalpha_dT

    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_dqs_dtem_dens_liq_0D( temp, dens,             & ! [IN]
                                                dqsat_dT_liq, qsat_liq_ ) ! [OUT]
    call ATMOS_SATURATION_dqs_dtem_dens_ice_0D( temp, dens,             & ! [IN]
                                                dqsat_dT_ice, qsat_ice_ ) ! [OUT]
    call ATMOS_SATURATION_alpha   ( temp, alpha_    ) ! [IN], [OUT]
    call ATMOS_SATURATION_dalphadT( temp, dalpha_dT ) ! [IN], [OUT]

    dqsat_dT  = qsat_liq_ * dalpha_dT + dqsat_dT_liq * (        alpha_ ) &
              - qsat_ice_ * dalpha_dT + dqsat_dT_ice * ( 1.0_RP-alpha_ )

    if ( present(qsat)     ) qsat = qsat_liq_ * alpha_ + qsat_ice_ * ( 1.0_RP-alpha_ )
    if ( present(qsat_liq) ) qsat_liq = qsat_liq_
    if ( present(qsat_ice) ) qsat_ice = qsat_ice_
    if ( present(alpha   ) ) alpha    = alpha_

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dens_all_0D

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  subroutine ATMOS_SATURATION_dqs_dtem_dpre_liq_0D( &
       temp, pres, qdry,   &
       dqsat_dT, dqsat_dP, &
       qsat, psat          )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    implicit none

    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: pres
    real(RP), intent(in)  :: qdry

    real(RP), intent(out) :: dqsat_dT
    real(RP), intent(out) :: dqsat_dP

    real(RP), intent(out), optional :: qsat
    real(RP), intent(out), optional :: psat

    real(RP) :: LHV ! latent heat of vaporization [J/kg]
    real(RP) :: psat_
    real(RP) :: den1, den2 ! denominator

    !---------------------------------------------------------------------------

    call HYDROMETEOR_LHV( temp, LHV )
    call ATMOS_SATURATION_psat_liq( temp, psat_ ) ! [IN], [OUT]

    den1 = ( pres - (1.0_RP-EPSvap) * psat_ )**2
    den2 = den1 * Rvap * temp**2

    dqsat_dP = - EPSvap * psat_ / den1
    dqsat_dT =   EPSvap * psat_ / den2 * LHV * pres

    if ( present(qsat) ) &
         call ATMOS_SATURATION_psat2qsat_pres( psat_, pres, qdry, & ! [IN]
                                               qsat               ) ! [OUT]
    if ( present(psat) ) psat = psat_

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dpre_liq_0D

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_SATURATION_dqs_dtem_dpre_liq_1D( &
       KA, KS, KE, &
       temp, pres, qdry,  &
       dqsat_dT, dqsat_dP )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp   (KA)
    real(RP), intent(in)  :: pres   (KA)
    real(RP), intent(in)  :: qdry   (KA)

    real(RP), intent(out) :: dqsat_dT(KA)
    real(RP), intent(out) :: dqsat_dP(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_dqs_dtem_dpre_liq_0D( temp(k), pres(k), qdry(k), & ! [IN]
                                                   dqsat_dT(k), dqsat_dP(k)   ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dpre_liq_1D

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  subroutine ATMOS_SATURATION_dqs_dtem_dpre_liq_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, pres, qdry,  &
       dqsat_dT, dqsat_dP )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: pres   (KA,IA,JA)
    real(RP), intent(in)  :: qdry   (KA,IA,JA)

    real(RP), intent(out) :: dqsat_dT(KA,IA,JA)
    real(RP), intent(out) :: dqsat_dP(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_dqs_dtem_dpre_liq_0D( temp(k,i,j), pres(k,i,j), qdry(k,i,j), & ! [IN]
                                                   dqsat_dT(k,i,j), dqsat_dP(k,i,j)       ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dpre_liq_3D

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqs_dtem_dpre_ice_0D( &
       temp, pres, qdry,   &
       dqsat_dT, dqsat_dP, &
       qsat                )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHS => ATMOS_HYDROMETEOR_LHS
    implicit none

    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: pres
    real(RP), intent(in)  :: qdry

    real(RP), intent(out) :: dqsat_dT
    real(RP), intent(out) :: dqsat_dP

    real(RP), intent(out), optional :: qsat

    real(RP) :: LHS ! latent heat of sublimation  [J/kg]
    real(RP) :: psat
    real(RP) :: den1, den2 ! denominator

    !---------------------------------------------------------------------------

    call HYDROMETEOR_LHS( temp, LHS ) ! [IN], [OUT]
    call ATMOS_SATURATION_psat_ice( temp, psat ) ! [IN], [OUT]

    den1 = ( pres - (1.0_RP-EPSvap) * psat )**2
    den2 = den1 * Rvap * temp**2

    dqsat_dP = - EPSvap * psat / den1
    dqsat_dT =   EPSvap * psat / den2 * LHS * pres

    if ( present(qsat) ) &
         call ATMOS_SATURATION_psat2qsat_pres( psat, pres, qdry, & ! [IN]
                                               qsat              ) ! [OUT]
    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dpre_ice_0D

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_SATURATION_dqs_dtem_dpre_ice_1D( &
       KA, KS, KE, &
       temp, pres, qdry,  &
       dqsat_dT, dqsat_dP )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp   (KA)
    real(RP), intent(in)  :: pres   (KA)
    real(RP), intent(in)  :: qdry   (KA)

    real(RP), intent(out) :: dqsat_dT(KA)
    real(RP), intent(out) :: dqsat_dP(KA)

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_dqs_dtem_dpre_ice( temp(k), pres(k), qdry(k), & ! [IN]
                                                dqsat_dT(k), dqsat_dP(k)   ) ! [OUT]
    enddo

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dpre_ice_1D

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqs_dtem_dpre_ice_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, pres, qdry,  &
       dqsat_dT, dqsat_dP )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: pres   (KA,IA,JA)
    real(RP), intent(in)  :: qdry   (KA,IA,JA)

    real(RP), intent(out) :: dqsat_dT(KA,IA,JA)
    real(RP), intent(out) :: dqsat_dP(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_dqs_dtem_dpre_ice( temp(k,i,j), pres(k,i,j), qdry(k,i,j), & ! [IN]
                                                dqsat_dT(k,i,j), dqsat_dP(k,i,j)       ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dpre_ice_3D

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqs_dtem_dpre_all_0D( &
       temp, pres, qdry,         &
       dqsat_dT, dqsat_dP,       &
       qsat, qsat_liq, qsat_ice, &
       alpha                     )
    implicit none

    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: pres
    real(RP), intent(in)  :: qdry

    real(RP), intent(out) :: dqsat_dT
    real(RP), intent(out) :: dqsat_dP

    real(RP), intent(out), optional :: qsat
    real(RP), intent(out), optional :: qsat_liq
    real(RP), intent(out), optional :: qsat_ice
    real(RP), intent(out), optional :: alpha

    real(RP) :: qsat_liq_, qsat_ice_, alpha_
    real(RP) :: dqsat_dT_liq, dqsat_dT_ice
    real(RP) :: dqsat_dP_liq, dqsat_dP_ice
    real(RP) :: dalpha_dT

    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_dqs_dtem_dpre_liq_0D( temp, pres, qdry,           & ! [IN]
                                                dqsat_dT_liq, dqsat_dP_liq, & ! [OUT]
                                                qsat_liq_                   ) ! [OUT]
    call ATMOS_SATURATION_dqs_dtem_dpre_ice_0D( temp, pres, qdry,           & ! [IN]
                                                dqsat_dT_ice, dqsat_dP_ice, & ! [OUT]
                                                qsat_ice_                   ) ! [OUT]
    call ATMOS_SATURATION_alpha   ( temp, alpha_    ) ! [IN], [OUT]
    call ATMOS_SATURATION_dalphadT( temp, dalpha_dT ) ! [IN], [OUT]

    dqsat_dT  = qsat_liq_ * dalpha_dT + dqsat_dT_liq * (        alpha_ ) &
              - qsat_ice_ * dalpha_dT + dqsat_dT_ice * ( 1.0_RP-alpha_ )
    dqsat_dP  = dqsat_dP_liq * (        alpha_ ) &
              + dqsat_dP_ice * ( 1.0_RP-alpha_ )

    if ( present(qsat)     ) qsat = qsat_liq * alpha_ + qsat_ice * ( 1.0_RP-alpha_ )
    if ( present(qsat_liq) ) qsat_liq = qsat_liq_
    if ( present(qsat_ice) ) qsat_ice = qsat_ice_
    if ( present(alpha   ) ) alpha    = alpha_

    return
  end subroutine ATMOS_SATURATION_dqs_dtem_dpre_all_0D

  !-----------------------------------------------------------------------------
  !> calculation of dew point
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_tdew_liq_0D( &
       DENS, TEMP, QV, &
       Tdew,           &
       converged       )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_LHV
    real(RP), intent(in) :: DENS
    real(RP), intent(in) :: TEMP
    real(RP), intent(in) :: QV

    real(RP), intent(out) :: Tdew
    logical,  intent(out) :: converged

    real(RP), parameter :: A = 17.625_RP
    real(RP), parameter :: B = 243.04_RP
    real(RP), parameter :: C = 610.94_RP
    integer,  parameter :: itelim = 100
    real(RP), parameter :: criteria = 0.1_RP**(2+RP/2)

    real(RP) :: lhv
    real(RP) :: pvap, psat
    real(RP) :: dpsat_dT
    real(RP) :: dTdew

    integer :: ite
    !---------------------------------------------------------------------------

    pvap = DENS * QV * Rvap * TEMP

    if ( pvap < psat_min_liq ) then
       converged = .true.
       Tdew = UNDEF
       return
    end if

    ! first guess is calculated by Alduchov and Eskridge (1996)
    ! See Lawrence (2005) BAMS
    Tdew = B * log( pvap / C ) / ( A - log( pvap / C ) ) + TEM00
    converged = .false.
    do ite = 1, itelim

       call ATMOS_SATURATION_psat_liq( Tdew, psat ) ! [IN], [OUT]
       call ATMOS_HYDROMETEOR_LHV( Tdew, lhv )

       dpsat_dT = psat * lhv / ( Rvap * Tdew**2 )
       dTdew = ( psat - pvap ) / dpsat_dT
       if ( dTdew < criteria ) then
          converged = .true.
          exit
       end if

       Tdew = Tdew - dTdew
    end do
    if( .not. converged ) then
       LOG_WARN("ATMOS_SATURATION_tdew_liq_0D",*) DENS, TEMP, QV, pvap, Tdew, dTdew, dpsat_dT
    endif

    return
  end subroutine ATMOS_SATURATION_tdew_liq_0D

!OCL SERIAL
  subroutine ATMOS_SATURATION_tdew_liq_1D( &
       KA, KS, KE, &
       DENS, TEMP, QV, &
       Tdew,           &
       converged       )
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: TEMP(KA)
    real(RP), intent(in) :: QV  (KA)

    real(RP), intent(out) :: Tdew(KA)
    logical,  intent(out) :: converged

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_SATURATION_tdew_liq_0D( DENS(k), TEMP(k), QV(k), & ! [IN]
                                          Tdew(k), converged       ) ! [OUT]
       if ( .not. converged ) exit
    end do

  end subroutine ATMOS_SATURATION_tdew_liq_1D

  subroutine ATMOS_SATURATION_tdew_liq_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, TEMP, QV, &
       Tdew            )
    use scale_prc, only: &
       PRC_abort
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: TEMP(KA,IA,JA)
    real(RP), intent(in) :: QV  (KA,IA,JA)

    real(RP), intent(out) :: Tdew(KA,IA,JA)

    logical :: converged, error
    integer :: k, i, j
    !---------------------------------------------------------------------------

    error = .false.
    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(converged)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_tdew_liq_0D( DENS(k,i,j), TEMP(k,i,j), QV(k,i,j), & ! [IN]
                                          Tdew(k,i,j), converged               ) ! [OUT]
       if ( .not. converged ) then
          LOG_ERROR("ATMOS_SATURATION_tdew_liq_3D",*) 'not converged! ', k,i,j
          error = .true.
          exit
       end if
    end do
    end do
    end do

    if ( error ) call PRC_abort

  end subroutine ATMOS_SATURATION_tdew_liq_3D

  !-----------------------------------------------------------------------------
  !> calculate equivalent potential temperature
  !>  Bolton, D., 1980: The computation of equivalent potential temperature. Monthly Weather Rev., 108, 1046-1053.
  !> PT_E = PT exp( L QV / (CPdry T) f )
  !> f ~ 1.0784 ( 1 + 0.810 QV )
  !> Here T_L is temperature at the lifting condensation level and
  !> T_L ~ 55 + 2840 / ( CPdry/Rdry log(T) - log(P_v) - 4.805 )
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_pote_0D( &
       DENS, POTT, TEMP, QV, &
       POTE                  )
    use scale_const, only: &
       Rdry  => CONST_Rdry, &
       CPdry => CONST_CPdry
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_LHV
    real(RP), intent(in) :: DENS
    real(RP), intent(in) :: POTT
    real(RP), intent(in) :: TEMP
    real(RP), intent(in) :: QV

    real(RP), intent(out) :: POTE

    real(RP) :: TL !> temperature at the lifting condensation level
    real(RP) :: Pv !> vapor pressure
    real(RP) :: LHV

    Pv = DENS * QV * Rvap * TEMP
    TL = 55.0_RP + 2840.0_RP / ( CPdry / Rdry * log(TEMP) - log(Pv) - 4.805_RP )
    call ATMOS_HYDROMETEOR_LHV( TEMP, LHV ) ! [IN], [OUT]

    POTE = POTT * exp( LHV * QV / ( CPdry * TEMP ) &
                     * 1.0784_RP * ( 1.0_RP + 0.810_RP * QV ) )

    return
  end subroutine ATMOS_SATURATION_pote_0D

!OCL SERIAL
  subroutine ATMOS_SATURATION_pote_1D( &
       KA, KS, KE, &
       DENS, POTT, TEMP, QV, &
       POTE                  )
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: POTT(KA)
    real(RP), intent(in) :: TEMP(KA)
    real(RP), intent(in) :: QV  (KA)

    real(RP), intent(out) :: POTE(KA)

    integer :: k

    do k = KS, KE
       call ATMOS_SATURATION_pote_0D( DENS(k), POTT(k), TEMP(k), QV(k), & ! [IN]
                                      POTE(k)                           ) ! [OUT]
    end do

    return
  end subroutine ATMOS_SATURATION_pote_1D

  subroutine ATMOS_SATURATION_pote_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, POTT, TEMP, QV, &
       POTE                  )
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: POTT(KA,IA,JA)
    real(RP), intent(in) :: TEMP(KA,IA,JA)
    real(RP), intent(in) :: QV  (KA,IA,JA)

    real(RP), intent(out) :: POTE(KA,IA,JA)

    integer :: k, i, j

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_SATURATION_pote_0D( &
            DENS(k,i,j), POTT(k,i,j), TEMP(k,i,j), QV(k,i,j), & ! [IN]
            POTE(k,i,j)                                       ) ! [OUT]
    end do
    end do
    end do

    return
  end subroutine ATMOS_SATURATION_pote_3D

  !-----------------------------------------------------------------------------
  !> Iterative moist conversion for liquid water at constant density (volume)
  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_SATURATION_moist_conversion_dens_liq_0D( &
       DENS, Emoist0,              &
       TEMP, QV, QC, CPtot, CVtot, &
       converged                   )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CP_VAPOR, &
       CP_WATER, &
       CV_VAPOR, &
       CV_WATER, &
       LHV
    implicit none

    real(RP), intent(in)    :: DENS
    real(RP), intent(in)    :: Emoist0

    real(RP), intent(inout) :: TEMP
    real(RP), intent(inout) :: QV
    real(RP), intent(inout) :: QC
    real(RP), intent(inout) :: CPtot
    real(RP), intent(inout) :: CVtot

    logical,  intent(out)   :: converged

    ! working
    real(RP) :: QSUM
    real(RP) :: CVtot0
    real(RP) :: qsat
    real(RP) :: Emoist ! moist internal energy

    ! d(X)/dT
    real(RP) :: dqsatl_dT
    real(RP) :: dqc_dT
    real(RP) :: dCVtot_dT
    real(RP) :: dEmoist_dT
    real(RP) :: dtemp

    integer,  parameter :: itelim = 100
    real(RP), parameter :: dtemp_criteria = 0.1_RP**(2+RP/2)
    integer  :: ite
    !---------------------------------------------------------------------------

    QSUM = QV + QC

    call ATMOS_SATURATION_dens2qsat_liq( temp, DENS, & ! [IN]
                                         qsat        ) ! [OUT]

    if ( QSUM <= qsat ) then
       QV = QSUM
       CPtot = CPtot + QC * ( CP_VAPOR - CP_WATER )
       CVtot = CVtot + QC * ( CV_VAPOR - CV_WATER )
       QC = 0.0_RP
       TEMP = ( Emoist0 - LHV * QV ) / CVtot
       converged = .true.
       return
    end if

    CVtot0 = CVtot

    converged = .false.
    do ite = 1, itelim

       ! get qsat and dX/dT under the temp
       call ATMOS_SATURATION_dqs_dtem_dens_liq_0D( temp, dens,     & ! [IN]
                                                   dqsatl_dT, qsat ) ! [OUT]

       dqc_dT = - dqsatl_dT

       dCVtot_dT = dqsatl_dT * CV_VAPOR &
                 + dqc_dT    * CV_WATER

       dEmoist_dT = temp * dCVtot_dT + CVtot + dqsatl_dT * LHV

       ! diagnose quantities with the qsat
       qc = QSUM - qsat
       CVtot = CVtot0 + ( CV_VAPOR - CV_WATER ) * ( qsat - QV )
       Emoist = temp * CVtot + qsat * LHV


       ! update temp by the newtonian method
       dtemp = ( Emoist - Emoist0 ) / dEmoist_dT
       temp  = temp - dtemp

       if ( abs(dtemp) < dtemp_criteria ) then
          converged = .true.
          exit
       endif

       if( temp*0.0_RP /= 0.0_RP ) exit
    enddo

    CPtot = CPtot + ( CP_VAPOR - CP_WATER ) * ( qsat - QV )

    QV = qsat

    return
  end subroutine ATMOS_SATURATION_moist_conversion_dens_liq_0d

  !-----------------------------------------------------------------------------
  !> Iterative moist conversion (liquid/ice mixture) at constant density (volume)
  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_SATURATION_moist_conversion_dens_all_0d( &
       DENS, Emoist0,                  &
       TEMP, QV, QC, QI, CPtot, CVtot, &
       converged                       )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CP_VAPOR, &
       CP_WATER, &
       CP_ICE, &
       CV_VAPOR, &
       CV_WATER, &
       CV_ICE, &
       LHV, &
       LHF
    implicit none

    real(RP), intent(in)    :: DENS
    real(RP), intent(in)    :: Emoist0

    real(RP), intent(inout) :: TEMP
    real(RP), intent(inout) :: QV
    real(RP), intent(inout) :: QC
    real(RP), intent(inout) :: QI
    real(RP), intent(inout) :: CPtot
    real(RP), intent(inout) :: CVtot

    logical,  intent(out)   :: converged

    ! working
    real(RP) :: QSUM
    real(RP) :: TEMP0
    real(RP) :: QV0, QC0, QI0
    real(RP) :: CVtot0
    real(RP) :: alpha
    real(RP) :: qsat
    real(RP) :: Emoist ! moist internal energy

    ! d(X)/dT
    real(RP) :: dalpha_dT
    real(RP) :: dqsat_dT
    real(RP) :: dqc_dT, dqi_dT
    real(RP) :: dCVtot_dT
    real(RP) :: dEmoist_dT
    real(RP) :: dtemp

    integer,  parameter :: itelim = 100
    real(RP), parameter :: dtemp_criteria = 0.1_RP**(2+RP/2)
    integer  :: ite
    !---------------------------------------------------------------------------

    TEMP0  = TEMP
    CVtot0 = CVtot
    QV0    = QV
    QC0    = QC
    QI0    = QI
    QSUM   = QV + QC + QI

    call ATMOS_SATURATION_dens2qsat_all( TEMP, DENS, & ! [IN]
                                         qsat        ) ! [OUT]

    if ( QSUM <= qsat ) then
       QV = QSUM
       CPtot = CPtot + QC * ( CP_VAPOR - CP_WATER ) + QI * ( CP_VAPOR - CP_ICE )
       CVtot = CVtot + QC * ( CV_VAPOR - CV_WATER ) + QI * ( CV_VAPOR - CV_ICE )
       QC = 0.0_RP
       QI = 0.0_RP
       TEMP = ( Emoist0 - LHV * QV ) / CVtot
       converged = .true.
       return
    end if

    converged = .false.
    do ite = 1, itelim

       ! dX/dT
       call ATMOS_SATURATION_dqs_dtem_dens_all_0D( temp, dens,            & ! [IN]
                                                   dqsat_dT,              & ! [OUT]
                                                   qsat=qsat, alpha=alpha ) ! [OUT]
       call ATMOS_SATURATION_dalphadT( temp, dalpha_dT ) ! [IN], [OUT]

       dqc_dT =  ( QSUM - qv ) * dalpha_dT - dqsat_dT * (        alpha )
       dqi_dT = -( QSUM - qv ) * dalpha_dT - dqsat_dT * ( 1.0_RP-alpha )

       dCVtot_dT = dqsat_dT * CV_VAPOR &
                 + dqc_dT   * CV_WATER &
                 + dqi_dT   * CV_ICE

       dEmoist_dT = temp * dCVtot_dT + CVtot + dqsat_dT * LHV - dqi_dT * LHF

       ! Saturation
       qv = qsat
       qc = ( QSUM - qsat ) * (          alpha )
       qi = ( QSUM - qsat ) * ( 1.0_RP - alpha )

       CVtot = CVtot0 &
             + CV_VAPOR * ( qv - QV0 ) &
             + CV_WATER * ( qc - QC0 ) &
             + CV_ICE   * ( qi - QI0 )

       Emoist = temp * CVtot + qv * LHV - qi * LHF

       ! update temp by the newtonian method
       dtemp = ( Emoist - Emoist0 ) / dEmoist_dT
       temp  = temp - dtemp

       if ( abs(dtemp) < dtemp_criteria ) then
          converged = .true.
          exit
       endif

       if( temp*0.0_RP /= 0.0_RP ) exit
    enddo

    CPtot = CPtot &
          + CP_VAPOR * ( qv - QV0 ) &
          + CP_WATER * ( qc - QC0 ) &
          + CP_ICE   * ( qi - QI0 )

    return
  end subroutine ATMOS_SATURATION_moist_conversion_dens_all_0d

  !-----------------------------------------------------------------------------
  !> Iterative moist conversion for liquid water at constant pressure
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_moist_conversion_pres_liq_0D( &
       PRES, Entr, Qdry, &
       QV, QC,           &
       Rtot, CPtot,      &
       TEMP,             &
       converged         )
    use scale_const, only: &
       EPS   => CONST_EPS, &
       LHV0  => CONST_LHV0
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_entr, &
       ATMOS_HYDROMETEOR_entr2temp, &
       CP_VAPOR, &
       CP_WATER

    real(RP), intent(in)    :: PRES
    real(RP), intent(in)    :: ENTR
    real(RP), intent(in)    :: Qdry

    real(RP), intent(inout) :: QV
    real(RP), intent(inout) :: QC
    real(RP), intent(inout) :: CPtot
    real(RP), intent(inout) :: Rtot

    real(RP), intent(out)   :: TEMP
    logical,  intent(out)   :: converged

    real(RP), parameter :: TEMMIN = 0.1_RP
    real(RP), parameter :: criteria = 1.E-8_RP
    integer,  parameter :: itelim   = 100
    integer             :: ite

    real(RP) :: Qsum
    real(RP) :: qsat, psat

    real(RP) :: TEMP_ite
    real(RP) :: QV_ite
    real(RP) :: ENTR_ite
    real(RP) :: Rtot_ite
    real(RP) :: CPtot_ite

    real(RP) :: dqsat_dT, dqsat_dP

    real(RP) :: TEMP_prev
    real(RP) :: dENTR_dT
    !---------------------------------------------------------------------------

    Qsum = QV + QC

    ! all liquid evaporates
    Rtot  = Rtot  + QC * Rvap
    CPtot = CPtot + QC * ( CP_VAPOR - CP_WATER )
    QV = Qsum
    QC = 0.0_RP

    call ATMOS_HYDROMETEOR_entr2temp( Entr, PRES, Qsum, 0.0_RP, Qdry, & ! [IN]
                                      Rtot, CPtot,                    & ! [IN]
                                      TEMP                            ) ! [OUT]

    call ATMOS_SATURATION_pres2qsat_liq( TEMP, PRES, Qdry, & ! [IN]
                                         qsat              ) ! [OUT]
    if ( Qsum <= qsat ) then
       ! unsaturated

       converged = .true.

       return

    else
       ! saturated

       TEMP_ite = TEMP

       do ite = 1, itelim

          call ATMOS_SATURATION_dqs_dtem_dpre_liq_0d( TEMP_ite, PRES, Qdry, & ! [IN]
                                                      dqsat_dT, dqsat_dP,   & ! [OUT]
                                                      qsat=qsat, psat=psat  ) ! [OUT]

          QV_ite = min( Qsum, qsat )

          Rtot_ite  = Rtot  - ( Qsum - QV_ite ) * Rvap
          CPtot_ite = CPtot - ( Qsum - QV_ite ) * ( CP_VAPOR - CP_WATER )

          dENTR_dT = CPtot_ite / TEMP_ite &
                   + ( ( CP_VAPOR - CP_WATER ) * log( TEMP_ite / TEM00 ) &
                     - Rvap * log( psat/ PSAT0 ) &
                     + LHV0 / TEM00 &
                     ) * dqsat_dT

          call ATMOS_HYDROMETEOR_entr( TEMP_ite, PRES,       & ! [IN]
                                       QV_ite, 0.0_RP, Qdry, & ! [IN] ! QI = 0
                                       Rtot_ite, CPtot_ite,  & ! [IN]
                                       ENTR_ite              ) ! [OUT]

          TEMP_prev = TEMP_ite
          TEMP_ite  = TEMP_ite - ( ENTR_ite - ENTR ) / max( dENTR_dT, EPS )
          TEMP_ite  = max( TEMP_ite, TEMMIN )

          if( abs(TEMP_ite-TEMP_prev) < criteria ) then
             converged = .true.
             exit
          end if

       enddo

    end if

    QV = QV_ite
    QC = Qsum - QV
    CPtot = CPtot_ite
    Rtot = Rtot_ite
    TEMP = TEMP_ite

    return
  end subroutine ATMOS_SATURATION_moist_conversion_pres_liq_0D

end module scale_atmos_saturation
