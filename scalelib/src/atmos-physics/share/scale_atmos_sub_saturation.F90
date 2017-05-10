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
     module procedure ATMOS_SATURATION_pres2qsat_all_3D_k
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
  real(RP), public :: CPovR_liq
  real(RP), public :: CPovR_ice
  real(RP), public :: CVovR_liq
  real(RP), public :: CVovR_ice
  real(RP), public :: LovR_liq
  real(RP), public :: LovR_ice

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
    implicit none

    NAMELIST / PARAM_ATMOS_SATURATION / &
       ATMOS_SATURATION_ULIMIT_TEMP, &
       ATMOS_SATURATION_LLIMIT_TEMP

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SATURATION] / Categ[ATMOS SHARE] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_SATURATION,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_SATURATION. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_SATURATION)

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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F7.2,A,F7.2)') '*** Temperature range for liquid/ice mixture : ', &
                                                      ATMOS_SATURATION_LLIMIT_TEMP, ' - ', &
                                                      ATMOS_SATURATION_ULIMIT_TEMP

    return
  end subroutine ATMOS_SATURATION_setup

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (0D)
  subroutine ATMOS_SATURATION_alpha_0D( &
       alpha, &
       temp   )
    implicit none

    real(RP), intent(out) :: alpha !< liquid/ice separation factor (0-1)
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

    real(RP), intent(out) :: alpha(KA) !< liquid/ice separation factor (0-1)
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

    real(RP), intent(out) :: alpha(KA,IA,JA) !< liquid/ice separation factor (0-1)
    real(RP), intent(in)  :: temp (KA,IA,JA) !< temperature [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JSB,JEB,ISB,IEB,KS,KE,alpha,temp,ATMOS_SATURATION_LLIMIT_TEMP,ATMOS_SATURATION_ULIMIT_TEMP)  
    do j = JSB, JEB
    do i = ISB, IEB
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

    call ATMOS_SATURATION_alpha   ( alpha, temp )
    call ATMOS_SATURATION_psat_liq( psatl, temp )
    call ATMOS_SATURATION_psat_ice( psati, temp )

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

    real(RP) :: alpha(KA), psatl(KA), psati(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_alpha   ( alpha(:), temp(:) )
    call ATMOS_SATURATION_psat_liq( psatl(:), temp(:) )
    call ATMOS_SATURATION_psat_ice( psati(:), temp(:) )

    do k = KS, KE
       psat(k) = psatl(k) * (          alpha(k) ) &
               + psati(k) * ( 1.0_RP - alpha(k) )
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

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,alpha,psatl,psati) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
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
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (0D)
  subroutine ATMOS_SATURATION_psat_liq_0D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp !< temperature               [K]
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( temp * RTEM00 )**CPovR_liq             &
                 * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp ) )

    return
  end subroutine ATMOS_SATURATION_psat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
  subroutine ATMOS_SATURATION_psat_liq_1D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA) !< temperature               [K]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       psat(k) = PSAT0 * ( temp(k) * RTEM00 )**CPovR_liq             &
                       * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k) ) )
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_liq_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (3D)
  subroutine ATMOS_SATURATION_psat_liq_3D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA,IA,JA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature               [K]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       psat(k,i,j) = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq             &
                           * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(k,i,j) ) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_liq_3D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (0D)
  subroutine ATMOS_SATURATION_psat_ice_0D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp !< temperature               [K]
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( temp * RTEM00 )**CPovR_ice             &
                 * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp ) )

    return
  end subroutine ATMOS_SATURATION_psat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
  subroutine ATMOS_SATURATION_psat_ice_1D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA) !< temperature               [K]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       psat(k) = PSAT0 * ( temp(k) * RTEM00 )**CPovR_ice             &
                       * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(k) ) )
    enddo

    return
  end subroutine ATMOS_SATURATION_psat_ice_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (3D)
  subroutine ATMOS_SATURATION_psat_ice_3D( &
       psat, &
       temp  )
    implicit none

    real(RP), intent(out) :: psat(KA,IA,JA) !< saturation vapor pressure [Pa]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature               [K]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       psat(k,i,j) = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice             &
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

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_all( psat, temp )

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

    real(RP) :: psat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_all( psat(:), temp(:) )

    do k = KS, KE
       qsat(k) = EPSvap * psat(k) / ( pres(k) - ( 1.0_RP-EPSvap ) * psat(k) )
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

    real(RP) :: psat

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       call ATMOS_SATURATION_psat_all( psat, temp(i,j) )

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

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,alpha,psatl,psati,psat) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
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
  !> calc temp & pres -> saturation vapor mass (liquid/ice mixture,3D)
  subroutine ATMOS_SATURATION_pres2qsat_all_3D_k( &
       qsat, &
       temp, &
       pres, &
       knum  )
    implicit none

    integer,  intent(in)  :: knum
    real(RP), intent(out) :: qsat(knum,IA,JA) !< saturation vapor mass [kg/kg]
    real(RP), intent(in)  :: temp(knum,IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: pres(knum,IA,JA) !< pressure              [Pa]

    real(RP) :: psat

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,psat) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = 1, knum
       call ATMOS_SATURATION_psat_all( psat, temp(k,i,j) )

       qsat(k,i,j) = EPSvap * psat / ( pres(k,i,j) - ( 1.0_RP-EPSvap ) * psat )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_pres2qsat_all_3D_k

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

    call ATMOS_SATURATION_psat_liq( psat, temp )

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

    real(RP) :: psat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_liq( psat(:), temp(:) )

    do k = KS, KE
       qsat(k) = EPSvap * psat(k) / ( pres(k) - ( 1.0_RP-EPSvap ) * psat(k) )
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

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,psat) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       psat = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq             &
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

    call ATMOS_SATURATION_psat_ice( psat, temp )

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

    real(RP) :: psat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_ice( psat(:), temp(:) )

    do k = KS, KE
       qsat(k) = EPSvap * psat(k) / ( pres(k) - ( 1.0_RP-EPSvap ) * psat(k) )
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

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,psat) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       psat = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice             &
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

    real(RP) :: psat
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_all( psat, temp )

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

    real(RP) :: psat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_all( psat(:), temp(:) )

    do k = KS, KE
       qsat(k) = psat(k) / ( dens(k) * Rvap * temp(k) )
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

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k,alpha,psatl,psati,psat) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JSB,JEB,ISB,IEB,KS,KE,temp,ATMOS_SATURATION_LLIMIT_TEMP,RTEM00,CPovR_liq,LovR_liq,qsat,dens,LovR_ice,CPovR_ice,ATMOS_SATURATION_ULIMIT_TEMP) 
    do j = JSB, JEB
    do i = ISB, IEB
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

    call ATMOS_SATURATION_psat_liq( psat, temp )

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

    real(RP) :: psat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_liq( psat(:), temp(:) )

    do k = KS, KE
       qsat(k) = psat(k) / ( dens(k) * Rvap * temp(k) )
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

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k,psat) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JSB,JEB,ISB,IEB,KS,KE,temp,RTEM00,CPovR_liq,LovR_liq,qsat,dens)  
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       psat = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_liq             &
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

    call ATMOS_SATURATION_psat_ice( psat, temp )

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

    real(RP) :: psat(KA)

    integer  :: k
    !---------------------------------------------------------------------------

    call ATMOS_SATURATION_psat_ice( psat(:), temp(:) )

    do k = KS, KE
       qsat(k) = psat(k) / ( dens(k) * Rvap * temp(k) )
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

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k,psat) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JSB,JEB,ISB,IEB,KS,KE,temp,RTEM00,CPovR_ice,LovR_ice,qsat,dens) 
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       psat = PSAT0 * ( temp(k,i,j) * RTEM00 )**CPovR_ice             &
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

    integer  :: k
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

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,lim1,lim2) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
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
  subroutine ATMOS_SATURATION_dqsw_dtem_rho( &
       dqsdtem, &
       temp,    &
       dens     )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA) !< (d qsw/d T)_{rho}
    real(RP), intent(in)  :: temp   (KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: dens   (KA,IA,JA) !< temperature           [K]

    real(RP) :: LHV(KA,IA,JA) ! latent heat of vaporization [J/kg]
    real(RP) :: psat
    real(RP) :: TEM

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    call HYDROMETEOR_LHV( LHV(:,:,:), temp(:,:,:) )

    !$omp parallel do private(i,j,k,TEM,psat) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat = PSAT0 * ( TEM * RTEM00 )**CPovR_liq             &
                    * exp( LovR_liq * ( RTEM00 - 1.0_RP/TEM ) )

       dqsdtem(k,i,j) = psat / ( dens(k,i,j) * Rvap * temp(k,i,j) * temp(k,i,j) ) &
                             * ( LHV(k,i,j) / ( Rvap * temp(k,i,j) ) - 1.0_RP )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqsw_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  subroutine ATMOS_SATURATION_dqsi_dtem_rho( &
       dqsdtem, &
       temp,    &
       dens     )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHS => ATMOS_HYDROMETEOR_LHS
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: dens   (KA,IA,JA)

    real(RP) :: LHS(KA,IA,JA) ! latent heat of sublimation  [J/kg]
    real(RP) :: psat
    real(RP) :: TEM

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    call HYDROMETEOR_LHS( LHS(:,:,:), temp(:,:,:) )

    !$omp parallel do private(i,j,k,TEM,psat) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat = PSAT0 * ( TEM * RTEM00 )**CPovR_ice             &
                    * exp( LovR_ice * ( RTEM00 - 1.0_RP/TEM ) )

       dqsdtem(k,i,j) = psat / ( dens(k,i,j) * Rvap * temp(k,i,j) * temp(k,i,j) ) &
                             * ( LHS(k,i,j) / ( Rvap * temp(k,i,j) ) - 1.0_RP )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqsi_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  subroutine ATMOS_SATURATION_dqsw_dtem_dpre( &
       dqsdtem, &
       dqsdpre, &
       temp,    &
       pres     )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(out) :: dqsdpre(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: pres   (KA,IA,JA)

    real(RP) :: LHV(KA,IA,JA) ! latent heat of vaporization [J/kg]
    real(RP) :: psat
    real(RP) :: TEM
    real(RP) :: den1, den2 ! denominator

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    call HYDROMETEOR_LHV( LHV(:,:,:), temp(:,:,:) )

    !$omp parallel do private(i,j,k,TEM,psat,den1,den2) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat = PSAT0 * ( TEM * RTEM00 )**CPovR_liq             &
                    * exp( LovR_liq * ( RTEM00 - 1.0_RP/TEM ) )

       den1 = ( pres(k,i,j) - (1.0_RP-EPSvap) * psat ) &
            * ( pres(k,i,j) - (1.0_RP-EPSvap) * psat )
       den2 = den1 * Rvap * temp(k,i,j) * temp(k,i,j)

       dqsdpre(k,i,j) = - EPSvap * psat / den1
       dqsdtem(k,i,j) =   EPSvap * psat / den2 * LHV(k,i,j) * pres(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqsw_dtem_dpre

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqsi_dtem_dpre( &
       dqsdtem, &
       dqsdpre, &
       temp,    &
       pres     )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHS => ATMOS_HYDROMETEOR_LHS
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(out) :: dqsdpre(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: pres   (KA,IA,JA)

    real(RP) :: LHS(KA,IA,JA) ! latent heat of sublimation  [J/kg]
    real(RP) :: psat
    real(RP) :: TEM
    real(RP) :: den1, den2 ! denominator

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    call HYDROMETEOR_LHS( LHS(:,:,:), temp(:,:,:) )

    !$omp parallel do private(i,j,k,TEM,psat,den1,den2) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat = PSAT0 * ( TEM * RTEM00 )**CPovR_ice             &
                    * exp( LovR_ice * ( RTEM00 - 1.0_RP/TEM ) )

       den1 = ( pres(k,i,j) - (1.0_RP-EPSvap) * psat ) &
            * ( pres(k,i,j) - (1.0_RP-EPSvap) * psat )
       den2 = den1 * Rvap * temp(k,i,j) * temp(k,i,j)

       dqsdpre(k,i,j) = - EPSvap * psat / den1
       dqsdtem(k,i,j) =   EPSvap * psat / den2 * LHS(k,i,j) * pres(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_SATURATION_dqsi_dtem_dpre

end module scale_atmos_saturation
