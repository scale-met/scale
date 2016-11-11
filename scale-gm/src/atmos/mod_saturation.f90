!-------------------------------------------------------------------------------
!> Module saturation process
!!
!! @par Description
!!         This module is for saturation processes
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_saturation
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
  public :: SATURATION_setup

  public :: SATURATION_alpha

  public :: SATURATION_psat_all
  public :: SATURATION_psat_liq
  public :: SATURATION_psat_ice

  interface SATURATION_alpha
     module procedure SATURATION_alpha_0D
     module procedure SATURATION_alpha_3D
  end interface SATURATION_alpha

  interface SATURATION_psat_all
     module procedure SATURATION_psat_all_0D
     module procedure SATURATION_psat_all_3D
  end interface SATURATION_psat_all
  interface SATURATION_psat_liq
     module procedure SATURATION_psat_liq_0D
     module procedure SATURATION_psat_liq_3D
  end interface SATURATION_psat_liq
  interface SATURATION_psat_ice
     module procedure SATURATION_psat_ice_0D
     module procedure SATURATION_psat_ice_3D
  end interface SATURATION_psat_ice

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public :: CPovR_liq
  real(RP), public :: CPovR_ice
  real(RP), public :: CVovR_liq
  real(RP), public :: CVovR_ice
  real(RP), public :: LovR_liq
  real(RP), public :: LovR_ice

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: TEM_MIN   = 10.0_RP !< minimum temperature [K]

  real(RP), private,      save :: SATURATION_ULIMIT_TEMP = 273.15_RP !< upper limit temperature
  real(RP), private,      save :: SATURATION_LLIMIT_TEMP = 233.15_RP !< lower limit temperature

  real(RP), private,      save :: RTEM00         !< inverse of TEM00
  real(RP), private,      save :: dalphadT_const !< d(alfa)/dt

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine SATURATION_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_THERMODYN_TYPE, &
       Rvap  => CONST_Rvap,  &
       CPvap => CONST_CPvap, &
       CVvap => CONST_CVvap, &
       CL    => CONST_CL,    &
       CI    => CONST_CI,    &
       LHV00 => CONST_LHV00, &
       LHS00 => CONST_LHS00, &
       LHV0  => CONST_LHV0,  &
       LHS0  => CONST_LHS0,  &
       TEM00 => CONST_TEM00
    implicit none

    NAMELIST / SATURATIONPARAM / &
       SATURATION_ULIMIT_TEMP, &
       SATURATION_LLIMIT_TEMP

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[saturation]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=SATURATIONPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** SATURATIONPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist SATURATIONPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=SATURATIONPARAM)

    RTEM00 = 1.0_RP / TEM00

    if ( CONST_THERMODYN_TYPE == 'EXACT' ) then

       CPovR_liq = ( CPvap - CL ) / Rvap
       CPovR_ice = ( CPvap - CI ) / Rvap
       CVovR_liq = ( CVvap - CL ) / Rvap
       CVovR_ice = ( CVvap - CI ) / Rvap

       LovR_liq  = LHV00 / Rvap
       LovR_ice  = LHS00 / Rvap

    elseif(      CONST_THERMODYN_TYPE == 'SIMPLE'  &
            .OR. CONST_THERMODYN_TYPE == 'SIMPLE2' ) then

       CPovR_liq = 0.0_RP
       CPovR_ice = 0.0_RP
       CVovR_liq = 0.0_RP
       CVovR_ice = 0.0_RP

       LovR_liq  = LHV0 / Rvap
       LovR_ice  = LHS0 / Rvap

    endif

    dalphadT_const = 1.0_RP / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F7.2,A,F7.2)') '*** Temperature range for ice : ', &
                                           SATURATION_LLIMIT_TEMP, ' - ', &
                                           SATURATION_ULIMIT_TEMP

    return
  end subroutine SATURATION_setup

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (0D)
  subroutine SATURATION_alpha_0D( &
       temp,  &
       alpha  )
    implicit none

    real(RP), intent(in)  :: temp  !< temperature [K]
    real(RP), intent(out) :: alpha !< liquid/ice separation factor [0-1]
    !---------------------------------------------------------------------------

    alpha = ( temp                   - SATURATION_LLIMIT_TEMP ) &
          / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
    alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

    return
  end subroutine SATURATION_alpha_0D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (3D)
  subroutine SATURATION_alpha_3D( &
       ijdim, &
       kdim,  &
       ldim,  &
       temp,  &
       alpha  )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: temp (ijdim,kdim,ldim) !< temperature [K]
    real(RP), intent(out) :: alpha(ijdim,kdim,ldim) !< liquid/ice separation factor [0-1]

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       alpha(ij,k,l) = ( temp(ij,k,l)           - SATURATION_LLIMIT_TEMP ) &
                     / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
       alpha(ij,k,l) = min( max( alpha(ij,k,l), 0.0_RP ), 1.0_RP )
    enddo
    enddo
    enddo

    return
  end subroutine SATURATION_alpha_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (0D)
  subroutine SATURATION_psat_all_0D( &
       temp,  &
       psat   )
    use scale_const, only: &
       PSAT0 => CONST_PSAT0
    implicit none

    real(RP), intent(in)  :: temp !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]

    real(RP) :: alpha, psatl, psati
    !---------------------------------------------------------------------------

    alpha = ( temp                   - SATURATION_LLIMIT_TEMP ) &
          / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
    alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

    psatl = PSAT0 * ( temp * RTEM00 )**CPovR_liq     &
          * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp ) )

    psati = PSAT0 * ( temp * RTEM00 )**CPovR_ice     &
          * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp ) )

    psat = psatl * (          alpha ) &
         + psati * ( 1.0_RP - alpha )

    return
  end subroutine SATURATION_psat_all_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (3D)
  subroutine SATURATION_psat_all_3D( &
       ijdim, &
       kdim,  &
       ldim,  &
       temp,  &
       psat   )
    use scale_const, only: &
       PSAT0 => CONST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: temp(ijdim,kdim,ldim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim,kdim,ldim) !< saturation vapor pressure [Pa]

    real(RP) :: alpha, psatl, psati

    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       alpha = ( temp(ij,k,l)           - SATURATION_LLIMIT_TEMP ) &
             / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       psatl = PSAT0 * ( temp(ij,k,l) * RTEM00 )**CPovR_liq             &
                     * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ij,k,l) ) )

       psati = PSAT0 * ( temp(ij,k,l) * RTEM00 )**CPovR_ice             &
                     * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ij,k,l) ) )

       psat(ij,k,l) = psatl * (          alpha ) &
                    + psati * ( 1.0_RP - alpha )
    enddo
    enddo
    enddo

    return
  end subroutine SATURATION_psat_all_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CL (0D)
  subroutine SATURATION_psat_liq_0D( &
       temp, &
       psat  )
    use scale_const, only: &
       PSAT0 => CONST_PSAT0
    implicit none

    real(RP), intent(in)  :: temp !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( temp * RTEM00 )**CPovR_liq     &
         * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp ) )

    return
  end subroutine SATURATION_psat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CL (3D)
  subroutine SATURATION_psat_liq_3D( &
       ijdim, &
       kdim,  &
       ldim,  &
       temp,  &
       psat   )
    use scale_const, only: &
       PSAT0 => CONST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: temp(ijdim,kdim,ldim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim,kdim,ldim) !< saturation vapor pressure [Pa]

    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       psat(ij,k,l) = PSAT0 * ( temp(ij,k,l) * RTEM00 )**CPovR_liq     &
                    * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ij,k,l) ) )
    enddo
    enddo
    enddo

    return
  end subroutine SATURATION_psat_liq_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CI (0D)
  subroutine SATURATION_psat_ice_0D( &
       temp, &
       psat  )
    use scale_const, only: &
       PSAT0 => CONST_PSAT0
    implicit none

    real(RP), intent(in)  :: temp !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( temp * RTEM00 )**CPovR_ice     &
         * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp ) )

    return
  end subroutine SATURATION_psat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CI (3D)
  subroutine SATURATION_psat_ice_3D( &
       ijdim, &
       kdim,  &
       ldim,  &
       temp,  &
       psat   )
    use scale_const, only: &
       PSAT0 => CONST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: temp(ijdim,kdim,ldim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim,kdim,ldim) !< saturation vapor pressure [Pa]

    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       psat(ij,k,l) = PSAT0 * ( temp(ij,k,l) * RTEM00 )**CPovR_ice     &
                    * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ij,k,l) ) )
    enddo
    enddo
    enddo

    return
  end subroutine SATURATION_psat_ice_3D

end module mod_saturation
