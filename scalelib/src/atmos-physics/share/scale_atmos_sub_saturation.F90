!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Saturation adjustment
!!
!! @par Description
!!          Saturation adjustment module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-24 (T.Seiki)   [new] Import from NICAM
!! @li      2012-02-10 (H.Yashiro) [mod] Reconstruction
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_saturation
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_const, only: &
     Rdry   => CONST_Rdry,   &
     CPdry  => CONST_CPdry,  &
     CVdry  => CONST_CVdry,  &
     Rvap   => CONST_Rvap,   &
     CPvap  => CONST_CPvap,  &
     CVvap  => CONST_CVvap,  &
     CL     => CONST_CL,     &
     CI     => CONST_CI,     &
     LHV00  => CONST_LH00,   &
     LHS00  => CONST_LHS00,  &
     LHV0   => CONST_LH0,    &
     LHS0   => CONST_LHS0,   &
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

  public :: ATMOS_SATURATION_pres2qsat_all
  public :: ATMOS_SATURATION_pres2qsat_liq
  public :: ATMOS_SATURATION_pres2qsat_ice

  public :: ATMOS_SATURATION_dens2qsat_all
  public :: ATMOS_SATURATION_dens2qsat_liq
  public :: ATMOS_SATURATION_dens2qsat_ice

  public :: ATMOS_SATURATION_dalphadT

  public :: ATMOS_SATURATION_dqsw_dtem_rho
  public :: ATMOS_SATURATION_dqsi_dtem_rho
  public :: ATMOS_SATURATION_dqsw_dtem_dpre
  public :: ATMOS_SATURATION_dqsi_dtem_dpre

  interface ATMOS_SATURATION_alpha
     module procedure ATMOS_SATURATION_alpha_0D
     module procedure ATMOS_SATURATION_alpha_1D
     module procedure ATMOS_SATURATION_alpha_3D
  end interface ATMOS_SATURATION_alpha

  interface ATMOS_SATURATION_psat_all
     module procedure ATMOS_SATURATION_psat_all_0D
     module procedure ATMOS_SATURATION_psat_all_1D
     module procedure ATMOS_SATURATION_psat_all_3D
  end interface ATMOS_SATURATION_psat_all
  interface ATMOS_SATURATION_psat_liq
     module procedure ATMOS_SATURATION_psat_liq_0D
     module procedure ATMOS_SATURATION_psat_liq_1D
     module procedure ATMOS_SATURATION_psat_liq_3D
  end interface ATMOS_SATURATION_psat_liq
  interface ATMOS_SATURATION_psat_ice
     module procedure ATMOS_SATURATION_psat_ice_0D
     module procedure ATMOS_SATURATION_psat_ice_1D
     module procedure ATMOS_SATURATION_psat_ice_3D
  end interface ATMOS_SATURATION_psat_ice

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

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: CPovR_liq
  real(RP), public, save :: CPovR_ice
  real(RP), public, save :: CVovR_liq
  real(RP), public, save :: CVovR_ice
  real(RP), public, save :: LovR_liq
  real(RP), public, save :: LovR_ice

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: TEM_MIN   = 10.0_RP !< minimum temperature [K]

  real(RP), private,      save :: ATMOS_SATURATION_ULIMIT_TEMP = 273.15_RP !< upper limit temperature
  real(RP), private,      save :: ATMOS_SATURATION_LLIMIT_TEMP = 233.15_RP !< lower limit temperature

  real(RP), private,      save :: RTEM00         !< inverse of TEM00
  real(RP), private,      save :: dalphadT_const !< d(alfa)/dt

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_SATURATION_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_SATURATION / &
       ATMOS_SATURATION_ULIMIT_TEMP, &
       ATMOS_SATURATION_LLIMIT_TEMP

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[SATURATION]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_SATURATION,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_SATURATION. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_SATURATION)

    RTEM00 = 1.0_RP / TEM00

    CPovR_liq = ( CPvap - CL ) / Rvap
    CPovR_ice = ( CPvap - CI ) / Rvap
    CVovR_liq = ( CVvap - CL ) / Rvap
    CVovR_ice = ( CVvap - CI ) / Rvap
    LovR_liq  = LHV00 / Rvap
    LovR_ice  = LHS00 / Rvap

    dalphadT_const = 1.0_RP / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )

    return
  end subroutine ATMOS_SATURATION_setup

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (0D)
  subroutine ATMOS_SATURATION_alpha_0D( &
       alpha, &
       temp   )
    implicit none

    real(RP), intent(out) :: alpha !< liquid/ice separation factor [0-1]
    real(RP), intent(in)  :: temp  !< temperature [K]
    !---------------------------------------------------------------------------

    alpha = ( temp                         - ATMOS_SATURATION_LLIMIT_TEMP ) &
          / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )

    alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

    return
  end subroutine ATMOS_SATURATION_alpha_0D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (1D)
  subroutine ATMOS_SATURATION_alpha_1D( &
       alpha, &
       temp   )
    implicit none

    real(RP), intent(out) :: alpha(KA) !< liquid/ice separation factor [0-1]
    real(RP), intent(in)  :: temp (KA) !< temperature [K]

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE

       alpha(k) = ( temp(k)                      - ATMOS_SATURATION_LLIMIT_TEMP ) &
                / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha(k) = min( max( alpha(k), 0.0_RP ), 1.0_RP )

    enddo

    return
  end subroutine ATMOS_SATURATION_alpha_1D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (3D)
  subroutine ATMOS_SATURATION_alpha_3D( &
       alpha, &
       temp   )
    implicit none

    real(RP), intent(out) :: alpha(KA,IA,JA) !< liquid/ice separation factor [0-1]
    real(RP), intent(in)  :: temp (KA,IA,JA) !< temperature [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       alpha(k,i,j) = ( temp(k,i,j)                  - ATMOS_SATURATION_LLIMIT_TEMP ) &
                    / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha(k,i,j) = min( max( alpha(k,i,j), 0.0_RP ), 1.0_RP )

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_alpha_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (0D)
  subroutine ATMOS_SATURATION_psat_all_0D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp !< temperature               [K]

    real(RP) :: alpha, psatl, psati
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_alpha_0D   ( alpha, temp )
    call ATMOS_SATURATION_psat_liq_0D( psatl, temp )
    call ATMOS_SATURATION_psat_ice_0D( psati, temp )

    psat = psatl * (          alpha ) &
         + psati * ( 1.0_RP - alpha )

    return
  end subroutine ATMOS_SATURATION_psat_all_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (1D)
  subroutine ATMOS_SATURATION_psat_all_1D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA) !< temperature               [K]

    real(RP) :: alpha, psatl, psati

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE

       alpha = ( temp(k)                      - ATMOS_SATURATION_LLIMIT_TEMP ) &
             / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       psatl = PSAT0 * ( temp(k) * RTEM00 )**CPovR_liq             &
                     * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k) ) )

       psati = PSAT0 * ( temp(k) * RTEM00 )**CPovR_ice             &
                     * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k) ) )

       psat(k) = psatl * (          alpha ) &
               + psati * ( 1.0_RP - alpha )
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_all_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (3D)
  subroutine ATMOS_SATURATION_psat_all_3D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA,IA,JA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature               [K]

    real(RP) :: alpha, psatl, psati

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,alpha,psatl,psati) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       alpha = ( temp(k,i,j)                  - ATMOS_SATURATION_LLIMIT_TEMP ) &
             / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       psatl = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq             &
                     * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       psati = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice             &
                     * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       psat(k,i,j) = psatl * (          alpha ) &
                   + psati * ( 1.0_RP - alpha )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_all_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CL (0D)
  subroutine ATMOS_SATURATION_psat_liq_0D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp !< temperature               [K]
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( temp * RTEM00 )**CPovR_liq     &
         * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp ) )

    return
  end subroutine ATMOS_SATURATION_psat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CL (1D)
  subroutine ATMOS_SATURATION_psat_liq_1D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA) !< temperature               [K]

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       psat(k) = PSAT0 * ( temp(k) * RTEM00 )**CPovR_liq     &
               * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k) ) )
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_liq_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CL (3D)
  subroutine ATMOS_SATURATION_psat_liq_3D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA,IA,JA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature               [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       psat(k,i,j) = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq     &
                   * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_liq_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CI (0D)
  subroutine ATMOS_SATURATION_psat_ice_0D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp !< temperature               [K]
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( temp * RTEM00 )**CPovR_ice     &
         * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp ) )

    return
  end subroutine ATMOS_SATURATION_psat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CI (1D)
  subroutine ATMOS_SATURATION_psat_ice_1D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA) !< temperature               [K]

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       psat(k) = PSAT0 * ( temp(k) * RTEM00 )**CPovR_ice     &
               * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k) ) )
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_ice_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation, based on CPV, CI (3D)
  subroutine ATMOS_SATURATION_psat_ice_3D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA,IA,JA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature               [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       psat(k,i,j) = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice     &
                   * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_ice_3D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,0D)
  subroutine ATMOS_SATURATION_pres2qsat_all_0D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp !< temperature           [K]
    real(RP), intent(in)  :: pres !< pressure              [Pa]

    real(RP) :: alpha, psatl, psati
    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_alpha_0D   ( alpha, temp )
    call ATMOS_SATURATION_psat_liq_0D( psatl, temp )
    call ATMOS_SATURATION_psat_ice_0D( psati, temp )

    psat = psatl * (          alpha ) &
         + psati * ( 1.0_RP - alpha )

    qsat = EPSvap * psat / ( pres - ( 1.0_RP-EPSvap ) * psat )

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_0D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,1D)
  subroutine ATMOS_SATURATION_pres2qsat_all_1D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat(KA) !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp(KA) !< temperature           [K]
    real(RP), intent(in)  :: pres(KA) !< pressure              [Pa]

    real(RP) :: alpha, psatl, psati
    real(RP) :: psat

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE

       alpha = ( temp(k)                      - ATMOS_SATURATION_LLIMIT_TEMP ) &
             / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       psatl = PSAT0 * ( temp(k) * RTEM00 )**CPovR_liq             &
                     * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k) ) )

       psati = PSAT0 * ( temp(k) * RTEM00 )**CPovR_ice             &
                     * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k) ) )

       psat = psatl * (          alpha ) &
            + psati * ( 1.0_RP - alpha )

       qsat(k) = EPSvap * psat / ( pres(k) - ( 1.0_RP-EPSvap ) * psat )

    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_1D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,2D)
  subroutine ATMOS_SATURATION_pres2qsat_all_2D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat(IA,JA) !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp(IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: pres(IA,JA) !< pressure              [Pa]

    real(RP) :: alpha, psatl, psati
    real(RP) :: psat

    integer :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,alpha,psatl,psati,psat) OMP_SCHEDULE_
    do j = JS, JE
    do i = IS, IE

       alpha = ( temp(i,j)                    - ATMOS_SATURATION_LLIMIT_TEMP ) &
             / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       psatl = PSAT0 * ( temp(i,j) * RTEM00 )**CPovR_liq             &
                     * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(i,j) ) )

       psati = PSAT0 * ( temp(i,j) * RTEM00 )**CPovR_ice             &
                     * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(i,j) ) )

       psat = psatl * (          alpha ) &
            + psati * ( 1.0_RP - alpha )

       qsat(i,j) = EPSvap * psat / ( pres(i,j) - ( 1.0_RP-EPSvap ) * psat )

    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_2D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,3D)
  subroutine ATMOS_SATURATION_pres2qsat_all_3D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat(KA,IA,JA) !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: pres(KA,IA,JA) !< pressure              [Pa]

    real(RP) :: alpha, psatl, psati
    real(RP) :: psat

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,alpha,psatl,psati,psat) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       alpha = ( temp(k,i,j)                  - ATMOS_SATURATION_LLIMIT_TEMP ) &
             / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       psatl = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq             &
                     * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       psati = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice             &
                     * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       psat = psatl * (          alpha ) &
            + psati * ( 1.0_RP - alpha )

       qsat(k,i,j) = EPSvap * psat / ( pres(k,i,j) - ( 1.0_RP-EPSvap ) * psat )

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_3D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid,0D)
  subroutine ATMOS_SATURATION_pres2qsat_liq_0D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp !< temperature           [K]
    real(RP), intent(in)  :: pres !< pressure              [Pa]

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_liq_0D( psat, temp )

    qsat = EPSvap * psat / ( pres - ( 1.0_RP-EPSvap ) * psat )

    return
  end subroutine ATMOS_SATURATION_pres2qsat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid,1D)
  subroutine ATMOS_SATURATION_pres2qsat_liq_1D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat(KA) !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp(KA) !< temperature           [K]
    real(RP), intent(in)  :: pres(KA) !< pressure              [Pa]

    real(RP) :: psat

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       psat = PSAT0 * ( temp(k) * RTEM00 )**CPovR_liq     &
            * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k) ) )

       qsat(k) = EPSvap * psat / ( pres(k) - ( 1.0_RP-EPSvap ) * psat )
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_liq_1D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (liquid,3D)
  subroutine ATMOS_SATURATION_pres2qsat_liq_3D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat(KA,IA,JA) !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: pres(KA,IA,JA) !< pressure              [Pa]

    real(RP) :: psat

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,psat) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       psat = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq     &
            * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       qsat(k,i,j) = EPSvap * psat / ( pres(k,i,j) - ( 1.0_RP-EPSvap ) * psat )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_liq_3D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (ice,0D)
  subroutine ATMOS_SATURATION_pres2qsat_ice_0D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp !< temperature           [K]
    real(RP), intent(in)  :: pres !< pressure              [Pa]

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_liq_0D( psat, temp )

    qsat = EPSvap * psat / ( pres - ( 1.0_RP-EPSvap ) * psat )

    return
  end subroutine ATMOS_SATURATION_pres2qsat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (ice,1D)
  subroutine ATMOS_SATURATION_pres2qsat_ice_1D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat(KA)
    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(in)  :: pres(KA)

    real(RP) :: psat

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       psat = PSAT0 * ( temp(k) * RTEM00 )**CPovR_ice     &
            * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k) ) )

       qsat(k) = EPSvap * psat / ( pres(k) - ( 1.0_RP-EPSvap ) * psat )
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_ice_1D

  !-----------------------------------------------------------------------------
  !> calc temp & pres -> saturation vapor mass (ice,3D)
  subroutine ATMOS_SATURATION_pres2qsat_ice_3D( &
       qsat, &
       temp, &
       pres  )
    implicit none

    real(RP), intent(out) :: qsat(KA,IA,JA)
    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: pres(KA,IA,JA)

    real(RP) :: psat

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,psat) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       psat = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice     &
            * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       qsat(k,i,j) = EPSvap * psat / ( pres(k,i,j) - ( 1.0_RP-EPSvap ) * psat )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_ice_3D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid/ice mixture,0D)
  subroutine ATMOS_SATURATION_dens2qsat_all_0D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat
    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP) :: alpha, psatl, psati
    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_alpha_0D   ( alpha, temp )
    call ATMOS_SATURATION_psat_liq_0D( psatl, temp )
    call ATMOS_SATURATION_psat_ice_0D( psati, temp )

    psat = psatl * (          alpha ) &
         + psati * ( 1.0_RP - alpha )

    qsat = psat / ( dens * Rvap * temp )

    return
  end subroutine ATMOS_SATURATION_dens2qsat_all_0D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid/ice mixture,1D)
  subroutine ATMOS_SATURATION_dens2qsat_all_1D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat(KA)
    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(in)  :: dens(KA)

    real(RP) :: alpha, psatl, psati
    real(RP) :: psat

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE

       alpha = ( temp(k)                      - ATMOS_SATURATION_LLIMIT_TEMP ) &
             / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       psatl = PSAT0 * ( temp(k) * RTEM00 )**CPovR_liq             &
                     * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k) ) )

       psati = PSAT0 * ( temp(k) * RTEM00 )**CPovR_ice             &
                     * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k) ) )

       psat = psatl * (          alpha ) &
            + psati * ( 1.0_RP - alpha )

       qsat(k) = psat / ( dens(k) * Rvap * temp(k) )

    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_all_1D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid/ice mixture,3D)
  subroutine ATMOS_SATURATION_dens2qsat_all_3D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat(KA,IA,JA)
    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: dens(KA,IA,JA)

    real(RP) :: alpha, psatl, psati
    real(RP) :: psat

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,alpha,psatl,psati,psat) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       alpha = ( temp(k,i,j)                  - ATMOS_SATURATION_LLIMIT_TEMP ) &
             / ( ATMOS_SATURATION_ULIMIT_TEMP - ATMOS_SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       psatl = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq             &
                     * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       psati = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice             &
                     * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       psat = psatl * (          alpha ) &
            + psati * ( 1.0_RP - alpha )

       qsat(k,i,j) = psat / ( dens(k,i,j) * Rvap * temp(k,i,j) )

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_all_3D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid,0D)
  subroutine ATMOS_SATURATION_dens2qsat_liq_0D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat
    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_liq_0D( psat, temp )

    qsat = psat / ( dens * Rvap * temp )

    return
  end subroutine ATMOS_SATURATION_dens2qsat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid,1D)
  subroutine ATMOS_SATURATION_dens2qsat_liq_1D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat(KA)
    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(in)  :: dens(KA)

    real(RP) :: psat

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       psat = PSAT0 * ( temp(k) * RTEM00 )**CPovR_liq     &
            * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k) ) )

       qsat(k) = psat / ( dens(k) * Rvap * temp(k) )
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_liq_1D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (liquid,3D)
  subroutine ATMOS_SATURATION_dens2qsat_liq_3D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat(KA,IA,JA)
    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: dens(KA,IA,JA)

    real(RP) :: psat

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,psat) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       psat = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq     &
            * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       qsat(k,i,j) = psat / ( dens(k,i,j) * Rvap * temp(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_liq_3D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (ice,0D)
  subroutine ATMOS_SATURATION_dens2qsat_ice_0D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat
    real(RP), intent(in)  :: temp
    real(RP), intent(in)  :: dens

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_ice_0D( psat, temp )

    qsat = psat / ( dens * Rvap * temp )

    return
  end subroutine ATMOS_SATURATION_dens2qsat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (ice,1D)
  subroutine ATMOS_SATURATION_dens2qsat_ice_1D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat(KA)
    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(in)  :: dens(KA)

    real(RP) :: psat

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       psat = PSAT0 * ( temp(k) * RTEM00 )**CPovR_ice     &
            * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k) ) )

       qsat(k) = psat / ( dens(k) * Rvap * temp(k) )
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_ice_1D

  !-----------------------------------------------------------------------------
  !> calc temp & dens -> saturation vapor mass (ice,3D)
  subroutine ATMOS_SATURATION_dens2qsat_ice_3D( &
       qsat, &
       temp, &
       dens  )
    implicit none

    real(RP), intent(out) :: qsat(KA,IA,JA)
    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: dens(KA,IA,JA)

    real(RP) :: psat

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,psat) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       psat = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice     &
            * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )

       qsat(k,i,j) = psat / ( dens(k,i,j) * Rvap * temp(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dens2qsat_ice_3D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(temp), 0D
  subroutine ATMOS_SATURATION_dalphadT_0D( &
       dalpha_dT, &
       temp       )
    implicit none

    real(RP), intent(out) :: dalpha_dT
    real(RP), intent(in)  :: temp

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
  subroutine ATMOS_SATURATION_dalphadT_1D( &
       dalpha_dT, &
       temp       )
    implicit none

    real(RP), intent(out) :: dalpha_dT(KA)
    real(RP), intent(in)  :: temp     (KA)

    real(RP) :: lim1, lim2

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE

       ! if Tup < temp(k), dalpha/dT = 0 (no slope)
       lim1 = 0.5_RP + sign( 0.5_RP, ATMOS_SATURATION_ULIMIT_TEMP - temp(k) )
       ! if Tdn > temp(k), dalpha/dT = 0 (no slope)
       lim2 = 0.5_RP + sign( 0.5_RP, temp(k) - ATMOS_SATURATION_LLIMIT_TEMP )

       dalpha_dT(k) = dalphadT_const * lim1 * lim2

    enddo

    return
  end subroutine ATMOS_SATURATION_dalphadT_1D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(temp), 3D
  subroutine ATMOS_SATURATION_dalphadT_3D( &
       dalpha_dT, &
       temp       )
    implicit none

    real(RP), intent(out) :: dalpha_dT(KA,IA,JA)
    real(RP), intent(in)  :: temp     (KA,IA,JA)

    real(RP) :: lim1, lim2

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,lim1,lim2) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       ! if Tup < temp(k,i,j), dalpha/dT = 0 (no slope)
       lim1 = 0.5_RP + sign( 0.5_RP, ATMOS_SATURATION_ULIMIT_TEMP - temp(k,i,j) )
       ! if Tdn > temp(k,i,j), dalpha/dT = 0 (no slope)
       lim2 = 0.5_RP + sign( 0.5_RP, temp(k,i,j) - ATMOS_SATURATION_LLIMIT_TEMP )

       dalpha_dT(k,i,j) = dalphadT_const * lim1 * lim2

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dalphadT_3D

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water
  subroutine ATMOS_SATURATION_dqsw_dtem_rho( dqsdtem, temp, dens )
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: dens   (KA,IA,JA)

    real(RP) :: psat(KA) ! saturation vapor pressure
    real(RP) :: lhv (KA) ! latent heat for condensation

    real(RP) :: RTEM00, CPovR, LovR, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    RTEM00   = 1.0_RP / TEM00
    CPovR = ( CPvap - CL ) / Rvap
    LovR = LHV00 / Rvap

    !$omp parallel do private(i,j,k,TEM,psat,lhv) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          TEM = max( temp(k,i,j), TEM_MIN )

          psat(k) = PSAT0                                  &
                   * ( TEM * RTEM00 )**CPovR            &
                   * exp( LovR * ( RTEM00 - 1.0_RP/TEM ) )
       enddo

       do k = KS, KE
          lhv(k)  = LHV0 + ( CPvap-CL ) * ( temp(k,i,j)-TEM00 )

          dqsdtem(k,i,j) = psat(k) / ( dens(k,i,j) * Rvap * temp(k,i,j) * temp(k,i,j) ) &
                         * ( lhv(k) / ( Rvap * temp(k,i,j) ) - 1.0_RP )
       enddo

    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqsw_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqsi_dtem_rho( dqsdtem, temp, dens )
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: dens   (KA,IA,JA)

    real(RP) :: psat(KA) ! saturation vapor pressure
    real(RP) :: lhv (KA) ! latent heat for condensation

    real(RP) :: RTEM00, CPovR, LovR, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    RTEM00   = 1.0_RP / TEM00
    CPovR = ( CPvap - CI ) / Rvap
    LovR = LHS00 / Rvap

    !$omp parallel do private(i,j,k,TEM,psat,lhv) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          TEM = max( temp(k,i,j), TEM_MIN )

          psat(k) = PSAT0                                 &
                  * ( TEM * RTEM00 )**CPovR            &
                  * exp( LovR * ( RTEM00 - 1.0_RP/TEM ) )
       enddo

       do k = KS, KE
          lhv(k) = LHS0 + ( CPvap-CI ) * ( temp(k,i,j)-TEM00 )

          dqsdtem(k,i,j) = psat(k) / ( dens(k,i,j) * Rvap * temp(k,i,j) * temp(k,i,j) ) &
                         * ( lhv(k) / ( Rvap * temp(k,i,j) ) - 1.0_RP )
       enddo

    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqsi_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqsw_dtem_dpre( dqsdtem, dqsdpre, temp, pres )
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(out) :: dqsdpre(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: pres   (KA,IA,JA)

    real(RP) :: psat(KA) ! saturation vapor pressure
    real(RP) :: lhv (KA) ! latent heat for condensation

    real(RP) :: den1(KA), den2(KA) ! denominator
    real(RP) :: RTEM00, CPovR, LovR, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    RTEM00   = 1.0_RP / TEM00
    CPovR = ( CPvap - CL ) / Rvap
    LovR = LHV00 / Rvap

    !$omp parallel do private(i,j,k,TEM,psat,den1,den2,lhv) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          TEM = max( temp(k,i,j), TEM_MIN )

          psat(k) = PSAT0                                 &
                  * ( TEM * RTEM00 )**CPovR            &
                  * exp( LovR * ( RTEM00 - 1.0_RP/TEM ) )
       enddo

       do k = KS, KE
          den1(k) = ( pres(k,i,j) - (1.0_RP-EPSvap) * psat(k) ) &
                  * ( pres(k,i,j) - (1.0_RP-EPSvap) * psat(k) )
          den2(k) = den1(k) * Rvap * temp(k,i,j) * temp(k,i,j)
          lhv (k) = LHV0 + ( CPvap-CL ) * ( temp(k,i,j)-TEM00 )
       enddo

       do k = KS, KE
          dqsdpre(k,i,j) = - EPSvap * psat(k) / den1(k)
          dqsdtem(k,i,j) =   EPSvap * psat(k) / den2(k) * lhv(k) * pres(k,i,j)
       enddo

    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqsw_dtem_dpre

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqsi_dtem_dpre( dqsdtem, dqsdpre, temp, pres )
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(out) :: dqsdpre(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: pres   (KA,IA,JA)

    real(RP) :: psat(KA) ! saturation vapor pressure
    real(RP) :: lhv (KA) ! latent heat for condensation

    real(RP) :: den1(KA), den2(KA) ! denominator
    real(RP) :: RTEM00, CPovR, LovR, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    RTEM00   = 1.0_RP / TEM00
    CPovR = ( CPvap - CI ) / Rvap
    LovR = LHS00 / Rvap

    !$omp parallel do private(i,j,k,TEM,psat,den1,den2,lhv) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          TEM = max( temp(k,i,j), TEM_MIN )

          psat(k) = PSAT0                                 &
                  * ( TEM * RTEM00 )**CPovR            &
                  * exp( LovR * ( RTEM00 - 1.0_RP/TEM ) )
       enddo

       do k = KS, KE
          den1(k) = ( pres(k,i,j) - (1.0_RP-EPSvap) * psat(k) ) &
                  * ( pres(k,i,j) - (1.0_RP-EPSvap) * psat(k) )
          den2(k) = den1(k) * Rvap * temp(k,i,j) * temp(k,i,j)
          lhv (k) = LHS0 + ( CPvap-CI ) * ( temp(k,i,j)-TEM00 )
       enddo

       do k = KS, KE
          dqsdpre(k,i,j) = - EPSvap * psat(k) / den1(k)
          dqsdtem(k,i,j) =   EPSvap * psat(k) / den2(k) * lhv(k) * pres(k,i,j)
       enddo

    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqsi_dtem_dpre

end module scale_atmos_saturation
!-------------------------------------------------------------------------------
