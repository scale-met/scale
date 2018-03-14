!-------------------------------------------------------------------------------
!> module Hydrometeor
!!
!! @par Description
!!          Hydrometeor module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-08-06 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_hydrometeor
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_HYDROMETEOR_setup
  public :: ATMOS_HYDROMETEOR_regist
  public :: ATMOS_HYDROMETEOR_LHV
  public :: ATMOS_HYDROMETEOR_LHS
  public :: ATMOS_HYDROMETEOR_LHF
  public :: ATMOS_HYDROMETEOR_entr
  public :: ATMOS_HYDROMETEOR_entr2temp
  public :: ATMOS_HYDROMETEOR_diagnose_number_concentration

  interface ATMOS_HYDROMETEOR_regist
     module procedure ATMOS_HYDROMETEOR_regist
     module procedure ATMOS_HYDROMETEOR_regist_obsolute
  end interface ATMOS_HYDROMETEOR_regist

  interface ATMOS_HYDROMETEOR_LHV
     module procedure ATMOS_HYDROMETEOR_LHV_0D
     module procedure ATMOS_HYDROMETEOR_LHV_1D
     module procedure ATMOS_HYDROMETEOR_LHV_2D
     module procedure ATMOS_HYDROMETEOR_LHV_3D
     module procedure ATMOS_HYDROMETEOR_LHV_1D_obsolute
     module procedure ATMOS_HYDROMETEOR_LHV_2D_obsolute
     module procedure ATMOS_HYDROMETEOR_LHV_3D_obsolute
  end interface ATMOS_HYDROMETEOR_LHV

  interface ATMOS_HYDROMETEOR_LHS
     module procedure ATMOS_HYDROMETEOR_LHS_0D
     module procedure ATMOS_HYDROMETEOR_LHS_1D
     module procedure ATMOS_HYDROMETEOR_LHS_2D
     module procedure ATMOS_HYDROMETEOR_LHS_3D
     module procedure ATMOS_HYDROMETEOR_LHS_3D_obsolute
  end interface ATMOS_HYDROMETEOR_LHS

  interface ATMOS_HYDROMETEOR_LHF
     module procedure ATMOS_HYDROMETEOR_LHF_0D
     module procedure ATMOS_HYDROMETEOR_LHF_1D
     module procedure ATMOS_HYDROMETEOR_LHF_2D
     module procedure ATMOS_HYDROMETEOR_LHF_3D
     module procedure ATMOS_HYDROMETEOR_LHF_3D_obsolute
  end interface ATMOS_HYDROMETEOR_LHF

  interface ATMOS_HYDROMETEOR_entr
     module procedure ATMOS_HYDROMETEOR_entr_0D
     module procedure ATMOS_HYDROMETEOR_entr_2D
     module procedure ATMOS_HYDROMETEOR_entr_3D
  end interface ATMOS_HYDROMETEOR_entr

  interface ATMOS_HYDROMETEOR_entr2temp
     module procedure ATMOS_HYDROMETEOR_entr2temp_0D
  end interface ATMOS_HYDROMETEOR_entr2temp

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public            :: ATMOS_HYDROMETEOR_ice_phase


  integer, public            :: I_QV = -1

  integer, public, parameter :: N_HYD = 6

  integer, public, parameter :: I_HC =  1 !< liquid water cloud
  integer, public, parameter :: I_HR =  2 !< liquid water rain
  integer, public, parameter :: I_HI =  3 !< ice water cloud
  integer, public, parameter :: I_HS =  4 !< snow
  integer, public, parameter :: I_HG =  5 !< graupel
  integer, public, parameter :: I_HH =  6 !< hail

  character(len=H_SHORT), public, parameter :: HYD_NAME(N_HYD) = &
       (/ "QC", "QR", "QI", "QS", "QG", "QH" /)
  character(len=H_MID), public, parameter :: HYD_DESC(N_HYD) = &
       (/ "cloud    ", "rain     ", "ice water", "snow     ", "graupel  ", "hail     " /)
  real(RP), public           :: HYD_DENS(N_HYD)

  integer, public            :: I_NC = -1
  integer, public            :: I_NR = -1
  integer, public            :: I_NI = -1
  integer, public            :: I_NS = -1
  integer, public            :: I_NG = -1
  integer, public            :: I_NH = -1

  integer, public            :: I_QC = -1
  integer, public            :: I_QR = -1
  integer, public            :: I_QI = -1
  integer, public            :: I_QS = -1
  integer, public            :: I_QG = -1
  integer, public            :: I_QH = -1

  ! hydrometeor (water + ice)
  integer, public            :: QHA =  0
  integer, public            :: QHS = -1
  integer, public            :: QHE = -2
  ! water
  integer, public            :: QLA =  0
  integer, public            :: QLS = -1
  integer, public            :: QLE = -2
  ! ice
  integer, public            :: QIA =  0
  integer, public            :: QIS = -1
  integer, public            :: QIE = -2

  real(RP), public           :: LHV       !< latent heat of vaporization for use [J/kg]
  real(RP), public           :: LHS       !< latent heat of sublimation  for use [J/kg]
  real(RP), public           :: LHF       !< latent heat of fusion       for use [J/kg]

  real(RP), public           :: CV_VAPOR !< CV for vapor [J/kg/K]
  real(RP), public           :: CP_VAPOR !< CP for vapor [J/kg/K]
  real(RP), public           :: CV_WATER !< CV for water [J/kg/K]
  real(RP), public           :: CP_WATER !< CP for water [J/kg/K]
  real(RP), public           :: CV_ICE   !< CV for ice   [J/kg/K]
  real(RP), public           :: CP_ICE   !< CP for ice   [J/kg/K]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: THERMODYN_EMASK = 1.0_RP

  logical,  private :: initialized = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_HYDROMETEOR_setup
    use scale_const, only: &
       CONST_setup, &
       CPvap          => CONST_CPvap, &
       CVvap          => CONST_CVvap, &
       CL             => CONST_CL,    &
       CI             => CONST_CI,    &
       LHV00          => CONST_LHV00, &
       LHS00          => CONST_LHS00, &
       LHF00          => CONST_LHF00, &
       LHV0           => CONST_LHV0,  &
       LHS0           => CONST_LHS0,  &
       LHF0           => CONST_LHF0,  &
       DWATR          => CONST_DWATR, &
       DICE           => CONST_DICE,  &
       THERMODYN_TYPE => CONST_THERMODYN_TYPE
    use scale_process, only: &
       PRC_abort
    implicit none
    !---------------------------------------------------------------------------

    if ( initialized ) return
    initialized = .true.

    call CONST_setup

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[HYDEROMETER] / Categ[ATMOS SHARE] / Origin[SCALElib]'

    if ( THERMODYN_TYPE == 'EXACT' ) then

       CV_VAPOR = CVvap
       CP_VAPOR = CPvap
       CV_WATER = CL
       CP_WATER = CV_WATER
       CV_ICE   = CI
       CP_ICE   = CV_ICE

       LHV = LHV00
       LHS = LHS00
       LHF = LHF00
       THERMODYN_EMASK = 1.0_RP

    elseif( THERMODYN_TYPE == 'SIMPLE' ) then

       CV_VAPOR = CVvap
       CP_VAPOR = CPvap
       CV_WATER = CVvap
       CP_WATER = CV_WATER
       CV_ICE   = CVvap
       CP_ICE   = CV_ICE

       LHV = LHV0
       LHS = LHS0
       LHF = LHF0
       THERMODYN_EMASK = 0.0_RP

    else
       write(*,*) 'xxx Not appropriate ATMOS_THERMODYN_ENERGY_TYPE. Check!', trim(THERMODYN_TYPE)
       call PRC_abort
    endif

    HYD_DENS(:) = (/ DWATR, & ! HC
                     DWATR, & ! HR
                     DICE,  & ! HI
                     DICE,  & ! HS
                     DICE,  & ! HG
                     DICE  /) ! HH

    return
  end subroutine ATMOS_HYDROMETEOR_setup

  !-----------------------------------------------------------------------------
  !> ATMOS_HYDROMETEOR_regist
  !! Regist tracer
  !<
  subroutine ATMOS_HYDROMETEOR_regist( &
       Q0,   &
       NL,   &
       NI,   &
       NAME, &
       DESC, &
       UNIT, &
       ADVC  )
    use scale_process, only: &
      PRC_abort
    use scale_tracer, only: &
      TRACER_regist
    use scale_const, only: &
      Rvap => CONST_Rvap
    implicit none

    integer,          intent(out) :: Q0
    integer,          intent(in)  :: NL             !< number of liquid water tracers
    integer,          intent(in)  :: NI             !< number of ice water tracers
    character(len=*), intent(in)  :: NAME(1+NL+NI)
    character(len=*), intent(in)  :: DESC(1+NL+NI)
    character(len=*), intent(in)  :: UNIT(1+NL+NI)

    logical,          intent(in), optional :: ADVC(1+NL+NI)

    real(RP) :: CV   (1+NL+NI)
    real(RP) :: CP   (1+NL+NI)
    real(RP) :: R    (1+NL+NI)
    logical  :: MASS (1+NL+NI)
    logical  :: ADVC_(1+NL+NI)

    integer  :: NQ
    integer  :: n
    !---------------------------------------------------------------------------

    if ( I_QV > 0 ) then
       write(*,*) 'xxx tracer for hydrometeor is already registerd'
       call PRC_abort
    endif

    NQ = 0

    ! vapor
    NQ = NQ + 1
    CV(NQ) = CV_VAPOR
    CP(NQ) = CP_VAPOR
    R (NQ) = Rvap

    ! liquid
    do n = 1, NL
       NQ = NQ + 1
       CV(NQ) = CV_WATER
       CP(NQ) = CP_WATER
       R (NQ) = 0.0_RP
    end do

    ! ice
    do n = 1, NI
       NQ = NQ + 1
       CV(NQ) = CV_ICE
       CP(NQ) = CP_ICE
       R (NQ) = 0.0_RP
    end do

    ! NQ = 1 + NL + NI,   vapor + liqid + ice

    if ( present(ADVC) ) then
       ADVC_(:) = ADVC(:)
    else
       ADVC_(:) = .true.
    endif

    do n = 1, NQ
       MASS(n) = .true.
    end do

    call TRACER_regist( Q0,         & ! [OUT]
                        NQ,         & ! [IN]
                        NAME,       & ! [IN]
                        DESC,       & ! [IN]
                        UNIT,       & ! [IN]
                        CV, CP, R,  & ! [IN], optional
                        ADVC_, MASS ) ! [IN], optional

    I_QV = Q0

    if ( NQ > 1 ) then
       QHS = I_QV + 1
       QHA = NL + NI
       QHE = QHS + QHA - 1
    endif

    if ( NL > 0 ) then
       QLS = I_QV + 1
       QLA = NL
       QLE = QLS + NL - 1
    endif

    if ( NI > 0 ) then
       QIS = QLE + 1
       QIA = NI
       QIE = QIS + NI - 1
    endif

    ATMOS_HYDROMETEOR_ice_phase = QIA > 0

    return
  end subroutine ATMOS_HYDROMETEOR_regist

  subroutine ATMOS_HYDROMETEOR_regist_obsolute( &
       Q0,   &
       NV,   &
       NL,   &
       NI,   &
       NAME, &
       DESC, &
       UNIT, &
       ADVC  )
    use scale_process, only: &
      PRC_MPIstop
    use scale_tracer, only: &
      TRACER_regist
    use scale_const, only: &
      Rvap => CONST_Rvap
    implicit none

    integer,          intent(out) :: Q0
    integer,          intent(in)  :: NV             !< number of vapor
    integer,          intent(in)  :: NL             !< number of liquid water tracers
    integer,          intent(in)  :: NI             !< number of ice water tracers
    character(len=*), intent(in)  :: NAME(NV+NL+NI)
    character(len=*), intent(in)  :: DESC(NV+NL+NI)
    character(len=*), intent(in)  :: UNIT(NV+NL+NI)

    logical,          intent(in), optional :: ADVC(NV+NL+NI)

    call ATMOS_HYDROMETEOR_regist( &
         Q0,               & ! [IN]
         NL, NI,           & ! [IN]
         NAME, DESC, UNIT, & ! [IN]
         ADVC              )

    return
  end subroutine ATMOS_HYDROMETEOR_regist_obsolute
  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_0D( &
       lhv, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHV0  => CONST_LHV0
    implicit none

    real(RP), intent(out) :: lhv  !< latent heat of vaporization [J/kg]
    real(RP), intent(in)  :: temp !< temperature                 [K]
    !---------------------------------------------------------------------------

    lhv = LHV0 + ( CP_VAPOR - CP_WATER ) * ( temp - TEM00 ) * THERMODYN_EMASK

    return
  end subroutine ATMOS_HYDROMETEOR_LHV_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_1D_obsolute( &
       lhv, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHV0  => CONST_LHV0
    implicit none

    real(RP), intent(out) :: lhv (KA) !< latent heat of vaporization [J/kg]
    real(RP), intent(in)  :: temp(KA) !< temperature                 [K]

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       lhv(k) = LHV0 + ( CP_VAPOR - CP_WATER ) * ( temp(k) - TEM00 ) * THERMODYN_EMASK
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHV_1D_obsolute

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_1D( &
       KA, KS, KE, &
       temp, &
       lhv   )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHV0  => CONST_LHV0
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA) !< temperature                 [K]

    real(RP), intent(out) :: lhv (KA) !< latent heat of vaporization [J/kg]

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       lhv(k) = LHV0 + ( CP_VAPOR - CP_WATER ) * ( temp(k) - TEM00 ) * THERMODYN_EMASK
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHV_1D
  
  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_2D_obsolute( &
       lhv, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHV0  => CONST_LHV0
    implicit none

    real(RP), intent(out) :: lhv (IA,JA) !< latent heat of vaporization [J/kg]
    real(RP), intent(in)  :: temp(IA,JA) !< temperature                 [K]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       lhv(i,j) = LHV0 + ( CP_VAPOR - CP_WATER ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHV_2D_obsolute

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_2D( &
       IA, IS, IE, JA, JS, JE, &
       temp, &
       lhv   )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHV0  => CONST_LHV0
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp(IA,JA) !< temperature                 [K]
    real(RP), intent(out) :: lhv (IA,JA) !< latent heat of vaporization [J/kg]

    integer :: i, j
    !---------------------------------------------------------------------------

    !omp paralell do
    do j = JS, JE
    do i = IS, IE
       lhv(i,j) = LHV0 + ( CP_VAPOR - CP_WATER ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHV_2D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_3D_obsolute( &
       lhv, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHV0  => CONST_LHV0
    implicit none

    real(RP), intent(out) :: lhv (KA,IA,JA) !< latent heat of vaporization [J/kg]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature                 [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       lhv(k,i,j) = LHV0 + ( CP_VAPOR - CP_WATER ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHV_3D_obsolute

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       lhv   )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHV0  => CONST_LHV0
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature                 [K]

    real(RP), intent(out) :: lhv (KA,IA,JA) !< latent heat of vaporization [J/kg]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       lhv(k,i,j) = LHV0 + ( CP_VAPOR - CP_WATER ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHV_3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHS_0D( &
       lhs, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHS0  => CONST_LHS0
    implicit none

    real(RP), intent(out) :: lhs   !< latent heat of sublimation [J/kg]
    real(RP), intent(in)  :: temp  !< temperature                [K]
    !---------------------------------------------------------------------------

    lhs = LHS0 + ( CP_VAPOR - CP_ICE ) * ( temp - TEM00 ) * THERMODYN_EMASK

    return
  end subroutine ATMOS_HYDROMETEOR_LHS_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHS_1D( &
       lhs, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHS0  => CONST_LHS0
    implicit none

    real(RP), intent(out) :: lhs (KA) !< latent heat of sublimation [J/kg]
    real(RP), intent(in)  :: temp(KA) !< temperature                [K]

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       lhs(k) = LHS0 + ( CP_VAPOR - CP_ICE ) * ( temp(k) - TEM00 ) * THERMODYN_EMASK
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHS_1D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHS_2D( &
       lhs, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHS0  => CONST_LHS0
    implicit none

    real(RP), intent(out) :: lhs (IA,JA) !< latent heat of sublimation [J/kg]
    real(RP), intent(in)  :: temp(IA,JA) !< temperature                [K]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       lhs(i,j) = LHS0 + ( CP_VAPOR - CP_ICE ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHS_2D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHS_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       lhs   )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHS0  => CONST_LHS0
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature                [K]
    real(RP), intent(out) :: lhs (KA,IA,JA) !< latent heat of sublimation [J/kg]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       lhs(k,i,j) = LHS0 + ( CP_VAPOR - CP_ICE ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHS_3D

  subroutine ATMOS_HYDROMETEOR_LHS_3D_obsolute( &
       lhs, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHS0  => CONST_LHS0
    implicit none

    real(RP), intent(out) :: lhs (KA,IA,JA) !< latent heat of sublimation [J/kg]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature                [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       lhs(k,i,j) = LHS0 + ( CP_VAPOR - CP_ICE ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHS_3D_obsolute

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHF_0D( &
       lhf, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHF0  => CONST_LHF0
    implicit none

    real(RP), intent(out) :: lhf  !< latent heat of fusion [J/kg]
    real(RP), intent(in)  :: temp !< temperature           [K]
    !---------------------------------------------------------------------------

    lhf = LHF0 + ( CP_WATER - CP_ICE ) * ( temp - TEM00 ) * THERMODYN_EMASK

    return
  end subroutine ATMOS_HYDROMETEOR_LHF_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHF_1D( &
       lhf, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHF0  => CONST_LHF0
    implicit none

    real(RP), intent(out) :: lhf (KA) !< latent heat of fusion [J/kg]
    real(RP), intent(in)  :: temp(KA) !< temperature           [K]

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       lhf(k) = LHF0 + ( CP_WATER - CP_ICE ) * ( temp(k) - TEM00 ) * THERMODYN_EMASK
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHF_1D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHF_2D( &
       lhf, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHF0  => CONST_LHF0
    implicit none

    real(RP), intent(out) :: lhf (IA,JA) !< latent heat of fusion [J/kg]
    real(RP), intent(in)  :: temp(IA,JA) !< temperature           [K]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       lhf(i,j) = LHF0 + ( CP_WATER - CP_ICE ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHF_2D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHF_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       temp, &
       lhf   )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHF0  => CONST_LHF0
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(out) :: lhf (KA,IA,JA) !< latent heat of fusion [J/kg]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       lhf(k,i,j) = LHF0 + ( CP_WATER - CP_ICE ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHF_3D

  subroutine ATMOS_HYDROMETEOR_LHF_3D_obsolute( &
       lhf, &
       temp )
    use scale_const, only: &
       TEM00 => CONST_TEM00, &
       LHF0  => CONST_LHF0
    implicit none

    real(RP), intent(out) :: lhf (KA,IA,JA) !< latent heat of fusion [J/kg]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature           [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       lhf(k,i,j) = LHF0 + ( CP_WATER - CP_ICE ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_LHF_3D_obsolute

  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (0D)
!OCL SERIAL
  subroutine ATMOS_HYDROMETEOR_entr_0D( &
       TEMP, PRES,   &
       QV, QI, Qdry, &
       Rtot, CPtot,  &
       entr          )
    use scale_const, only: &
       PRE00 => CONST_PRE00, &
       TEM00 => CONST_TEM00, &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       LHV0  => CONST_LHV0,  &
       LHF0  => CONST_LHF0,  &
       PSAT0 => CONST_PSAT0
    implicit none

    real(RP), intent(in) :: TEMP  !< temperature [K]
    real(RP), intent(in) :: PRES  !< pressure    [Pa]
    real(RP), intent(in) :: QV    !< water vapor  mass concentration [kg/kg]
    real(RP), intent(in) :: QI    !< ice water    mass concentration [kg/kg]
    real(RP), intent(in) :: Qdry  !< dry air      mass concentration [kg/kg]
    real(RP), intent(in) :: Rtot  !< gas constant
    real(RP), intent(in) :: CPtot !< specific heat

    real(RP), intent(out) :: entr !< entropy            [J/K]

    real(RP) :: pres_dry, pres_vap
    !---------------------------------------------------------------------------

    pres_dry = PRES * Qdry * Rdry / Rtot
    pres_vap = PRES * QV   * Rvap / Rtot

    entr = CPtot * log( TEMP / TEM00 ) &
         - Qdry * Rdry * log( pres_dry / PRE00 ) &
         - QV   * Rvap * log( pres_vap / PSAT0 ) &
         + ( QV * LHV0 - QI * LHF0 ) / TEM00

    return
  end subroutine ATMOS_HYDROMETEOR_entr_0D

  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (2D)
  subroutine ATMOS_HYDROMETEOR_entr_2D( &
       IA, IS, IE, JA, JS, JE, &
       TEMP, PRES,   &
       QV, QI, Qdry, &
       Rtot, CPtot,  &
       entr          )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: TEMP (IA,JA) !< temperature [K]
    real(RP), intent(in) :: PRES (IA,JA) !< pressure    [Pa]
    real(RP), intent(in) :: QV   (IA,JA) !< water vapor  mass concentration [kg/kg]
    real(RP), intent(in) :: QI   (IA,JA) !< ice water    mass concentration [kg/kg]
    real(RP), intent(in) :: Qdry (IA,JA) !< dry air      mass concentration [kg/kg]
    real(RP), intent(in) :: Rtot (IA,JA) !< gas constant
    real(RP), intent(in) :: CPtot(IA,JA) !< specific heat

    real(RP), intent(out) :: entr(IA,JA)    !< entropy            [J/K]

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROMETEOR_entr_0D( &
            TEMP(i,j), PRES(i,j),        & ! [IN]
            QV(i,j), QI(i,j), Qdry(i,j), & ! [IN]
            Rtot(i,j), CPtot(i,j),       & ! [IN]
            entr(i,j)                    ) ! [OUT]
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_entr_2D
  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (3D)
  subroutine ATMOS_HYDROMETEOR_entr_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       TEMP, PRES,   &
       QV, QI, Qdry, &
       Rtot, CPtot,  &
       entr          )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: TEMP (KA,IA,JA) !< temperature [K]
    real(RP), intent(in) :: PRES (KA,IA,JA) !< pressure    [Pa]
    real(RP), intent(in) :: QV   (KA,IA,JA) !< water vapor  mass concentration [kg/kg]
    real(RP), intent(in) :: QI   (KA,IA,JA) !< ice water    mass concentration [kg/kg]
    real(RP), intent(in) :: Qdry (KA,IA,JA) !< dry air      mass concentration [kg/kg]
    real(RP), intent(in) :: Rtot (KA,IA,JA) !< gas constant
    real(RP), intent(in) :: CPtot(KA,IA,JA) !< specific heat

    real(RP), intent(out) :: entr(KA,IA,JA) !< entropy            [J/K]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_HYDROMETEOR_entr_0D( &
            TEMP(k,i,j), PRES(k,i,j),          & ! [IN]
            QV(k,i,j), QI(k,i,j), Qdry(k,i,j), & ! [IN]
            Rtot(k,i,j), CPtot(k,i,j),         & ! [IN]
            entr(k,i,j)                        ) ! [OUT]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_entr_3D

  !-----------------------------------------------------------------------------
  !> calc entropy, pres, q -> temp (0D)
!OCL SERIAL
  subroutine ATMOS_HYDROMETEOR_entr2temp_0D( &
       entr, pres, &
       qv, qi, qdry, &
       Rtot, CPtot,    &
       temp            )
    use scale_const, only: &
       PRE00 => CONST_PRE00, &
       TEM00 => CONST_TEM00, &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       LHV0  => CONST_LHV0,  &
       LHF0  => CONST_LHF0,  &
       PSAT0 => CONST_PSAT0
    implicit none

    real(RP), intent(in) :: entr  !< entropy  [J/K]
    real(RP), intent(in) :: pres  !< pressure [Pa]
    real(RP), intent(in) :: qv    !< water vapor mass concentration [kg/kg]
    real(RP), intent(in) :: qi    !< ice water   mass concentration [kg/kg]
    real(RP), intent(in) :: qdry  !< dry air     mass concentration [kg/kg]
    real(RP), intent(in) :: Rtot  !< total gas constant
    real(RP), intent(in) :: CPtot !< total specific heat

    real(RP), intent(out) :: temp !< temperature [K]

    real(RP) :: pres_dry, pres_vap
    !---------------------------------------------------------------------------

    pres_dry = PRES * Qdry * Rdry / Rtot
    pres_vap = PRES * QV   * Rvap / Rtot

    TEMP = TEM00 &
         * exp( ( entr &
                + Qdry * Rdry * log( pres_dry / PRE00 ) &
                + QV   * Rvap * log( pres_vap / PSAT0 ) &
              - ( QV * LHV0 - QI * LHF0 ) / TEM00 ) / CPtot )
    return
  end subroutine ATMOS_HYDROMETEOR_entr2temp_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_diagnose_number_concentration( &
       QTRC )
    use scale_const, only: &
       PI => CONST_PI
    implicit none

    real(RP), intent(inout) :: QTRC(:,:,:,:)

    real(RP), parameter :: Dc   =  20.E-6_RP ! typical particle diameter for cloud  [m]
    real(RP), parameter :: Dr   = 200.E-6_RP ! typical particle diameter for rain   [m]
    real(RP), parameter :: Di   =  80.E-6_RP ! typical particle diameter for ice    [m]
    real(RP), parameter :: Ds   =  80.E-6_RP ! typical particle diameter for snow   [m]
    real(RP), parameter :: Dg   = 200.E-6_RP ! typical particle diameter for grapel [m]
    real(RP), parameter :: RHOw =  1000.0_RP ! typical density for warm particles   [kg/m3]
    real(RP), parameter :: RHOf =   100.0_RP ! typical density for frozen particles [kg/m3]
    real(RP), parameter :: RHOg =   400.0_RP ! typical density for grapel particles [kg/m3]
    real(RP), parameter :: b    =     3.0_RP ! assume spherical form

    real(RP) :: piov6
    !---------------------------------------------------------------------------

    piov6 = pi / 6.0_RP

    if ( I_NC > 0 ) QTRC(:,:,:,I_NC) = QTRC(:,:,:,I_QC) / ( (piov6*RHOw) * Dc**b )
    if ( I_NR > 0 ) QTRC(:,:,:,I_NR) = QTRC(:,:,:,I_QR) / ( (piov6*RHOw) * Dr**b )
    if ( I_NI > 0 ) QTRC(:,:,:,I_NI) = QTRC(:,:,:,I_QI) / ( (piov6*RHOf) * Di**b )
    if ( I_NS > 0 ) QTRC(:,:,:,I_NS) = QTRC(:,:,:,I_QS) / ( (piov6*RHOf) * Ds**b )
    if ( I_NG > 0 ) QTRC(:,:,:,I_NG) = QTRC(:,:,:,I_QG) / ( (piov6*RHOg) * Dg**b )

    return
  end subroutine ATMOS_HYDROMETEOR_diagnose_number_concentration

end module scale_atmos_hydrometeor
