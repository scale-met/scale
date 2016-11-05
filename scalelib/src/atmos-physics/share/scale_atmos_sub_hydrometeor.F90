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
module scale_atmos_hydrometeor
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
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
  public :: ATMOS_HYDROMETEOR_diagnose_number_concentration

  interface ATMOS_HYDROMETEOR_LHV
     module procedure ATMOS_HYDROMETEOR_LHV_0D
     module procedure ATMOS_HYDROMETEOR_LHV_1D
     module procedure ATMOS_HYDROMETEOR_LHV_2D
     module procedure ATMOS_HYDROMETEOR_LHV_3D
  end interface ATMOS_HYDROMETEOR_LHV

  interface ATMOS_HYDROMETEOR_LHS
     module procedure ATMOS_HYDROMETEOR_LHS_0D
     module procedure ATMOS_HYDROMETEOR_LHS_1D
     module procedure ATMOS_HYDROMETEOR_LHS_2D
     module procedure ATMOS_HYDROMETEOR_LHS_3D
  end interface ATMOS_HYDROMETEOR_LHS

  interface ATMOS_HYDROMETEOR_LHF
     module procedure ATMOS_HYDROMETEOR_LHF_0D
     module procedure ATMOS_HYDROMETEOR_LHF_1D
     module procedure ATMOS_HYDROMETEOR_LHF_2D
     module procedure ATMOS_HYDROMETEOR_LHF_3D
  end interface ATMOS_HYDROMETEOR_LHF

  interface ATMOS_HYDROMETEOR_entr
     module procedure ATMOS_HYDROMETEOR_entr_0D
     module procedure ATMOS_HYDROMETEOR_entr_2D
     module procedure ATMOS_HYDROMETEOR_entr_3D
  end interface ATMOS_HYDROMETEOR_entr

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public            :: I_QV = -1

  integer, public, parameter :: N_HYD = 6

  integer, public, parameter :: I_HC =  1 !< liquid water cloud
  integer, public, parameter :: I_HR =  2 !< liquid water rain
  integer, public, parameter :: I_HI =  3 !< ice water cloud
  integer, public, parameter :: I_HS =  4 !< snow
  integer, public, parameter :: I_HG =  5 !< graupel
  integer, public, parameter :: I_HH =  6 !< hail

  integer, public            :: I_QC = -1
  integer, public            :: I_QR = -1
  integer, public            :: I_QI = -1
  integer, public            :: I_QS = -1
  integer, public            :: I_QG = -1
  integer, public            :: I_QH = -1

  integer, public            :: I_NC = -1
  integer, public            :: I_NR = -1
  integer, public            :: I_NI = -1
  integer, public            :: I_NS = -1
  integer, public            :: I_NG = -1
  integer, public            :: I_NH = -1

  ! hydrometeor (water + ice)
  integer, public            :: QHS = -1
  integer, public            :: QHE = -2
  ! water
  integer, public            :: QLS = -1
  integer, public            :: QLE = -2
  ! ice
  integer, public            :: QIS = -1
  integer, public            :: QIE = -2

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: CV_VAPOR !< CV for vapor [J/kg/K]
  real(RP), private :: CP_VAPOR !< CP for vapor [J/kg/K]
  real(RP), private :: CV_WATER !< CV for water [J/kg/K]
  real(RP), private :: CP_WATER !< CP for water [J/kg/K]
  real(RP), private :: CV_ICE   !< CV for ice   [J/kg/K]
  real(RP), private :: CP_ICE   !< CP for ice   [J/kg/K]

  real(RP), private :: THERMODYN_EMASK = 1.0_RP

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_HYDROMETEOR_setup
    use scale_const, only: &
       CPdry          => CONST_CPdry,          &
       CVdry          => CONST_CVdry,          &
       CPvap          => CONST_CPvap,          &
       CVvap          => CONST_CVvap,          &
       CL             => CONST_CL,             &
       CI             => CONST_CI,             &
       LHV00          => CONST_LHV00,          &
       LHS00          => CONST_LHS00,          &
       LHF00          => CONST_LHF00,          &
       LHV0           => CONST_LHV0,           &
       LHS0           => CONST_LHS0,           &
       LHF0           => CONST_LHF0,           &
       THERMODYN_TYPE => CONST_THERMODYN_TYPE
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    !---------------------------------------------------------------------------

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
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_HYDROMETEOR_setup

  !-----------------------------------------------------------------------------
  !> Regist tracer
  subroutine ATMOS_HYDROMETEOR_regist( &
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

    real(RP) :: CV(NV+NL+NI)
    real(RP) :: CP(NV+NL+NI)
    real(RP) :: R(NV+NL+NI)
    logical  :: MASS(NV+NL+NI)
    logical  :: ADVC_(NV+NL+NI)

    integer  :: NQ
    integer  :: n
    !---------------------------------------------------------------------------

    if ( I_QV > 0 ) then
       write(*,*) 'xxx tracer for hydrometeor is already registerd'
       call PRC_MPIstop
    endif

    if ( NV /= 1 ) then
       write(*,*) 'xxx number of vapor must be 1 at this moment'
       call PRC_MPIstop
    endif

    NQ = 0

    do n = 1, NV
       NQ = NQ + 1
       CV(NQ) = CV_VAPOR
       CP(NQ) = CP_VAPOR
       R (NQ) = Rvap
    end do

    do n = 1, NL
       NQ = NQ + 1
       CV(NQ) = CV_WATER
       CP(NQ) = CP_WATER
       R (NQ) = 0.0_RP
    end do

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

    call TRACER_regist( Q0,                    & ! (out)
                        NQ, NAME, DESC, UNIT,  & ! (in)
                        CV, CP, R, ADVC_, MASS ) ! (in)

    I_QV = Q0

    if ( NQ > 1 ) then
       QHS = I_QV + 1
       QHE = QHS + NL + NI - 1
    endif

    if ( NL > 0 ) then
       QLS = I_QV + 1
       QLE = QLS + NL - 1
    endif

    if ( NI > 0 ) then
       QIS = QLE + 1
       QIE = QIS + NI - 1
    endif

    return
  end subroutine ATMOS_HYDROMETEOR_regist

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
  subroutine ATMOS_HYDROMETEOR_LHV_1D( &
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
  end subroutine ATMOS_HYDROMETEOR_LHV_1D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_2D( &
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
  end subroutine ATMOS_HYDROMETEOR_LHV_2D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_LHV_3D( &
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
  end subroutine ATMOS_HYDROMETEOR_LHS_3D

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
  end subroutine ATMOS_HYDROMETEOR_LHF_3D

  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (0D)
  subroutine ATMOS_HYDROMETEOR_entr_0D( &
       entr, &
       temp, &
       pres, &
       q,    &
       Rq    )
    use scale_const, only: &
       EPS   => CONST_EPS,   &
       PRE00 => CONST_PRE00, &
       TEM00 => CONST_TEM00, &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       Rvap  => CONST_Rvap,  &
       LHV0  => CONST_LHV0,  &
       LHF0  => CONST_LHF0,  &
       PSAT0 => CONST_PSAT0
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(out) :: entr     !< entropy            [J/K]
    real(RP), intent(in)  :: temp     !< temperature        [K]
    real(RP), intent(in)  :: pres     !< pressure           [Pa]
    real(RP), intent(in)  :: q(QA)    !< mass concentration [kg/kg]
    real(RP), intent(in)  :: Rq(QA)   !< gas constantt      [J/kg/K]

    real(RP) :: qdry, qliq, qice, Rtot
    real(RP) :: logT_T0, pres_dry, pres_vap

    integer  :: iqw
    !---------------------------------------------------------------------------

    logT_T0 = log( temp / TEM00 )

    qdry = 1.0_RP
    Rtot = 0.0_RP
    do iqw = 1, QA
       qdry = qdry - q(iqw)
       Rtot = Rtot + q(iqw) * Rq(iqw)
    enddo
    Rtot = Rtot + Rdry * qdry

    ! dry air + vapor
    pres_dry = max( pres * qdry * Rdry / Rtot, EPS )
    entr = qdry * CPdry * logT_T0 &
         - qdry * Rdry  * log( pres_dry / PRE00 )

    if ( I_QV > 0 ) then
       pres_vap = max( pres * q(I_QV) * Rvap / Rtot, EPS )
       entr = entr + q(I_QV) * CP_VAPOR * logT_T0 &
                   - q(I_QV) * Rvap  * log( pres_vap / PSAT0 ) &
                   + q(I_QV) * LHV0 / TEM00
    endif

    ! liquid water
    if ( QLS > 0 ) then
       do iqw = QLS, QLE
          entr = entr + q(iqw) * CP_WATER * logT_T0
       enddo
    endif

    ! ice
    if ( QIS > 0 ) then
       do iqw = QIS, QIE
          entr = entr + q(iqw) * CP_ICE * logT_T0 &
                      - q(iqw) * LHF0 / TEM00
       enddo
    endif

    return
  end subroutine ATMOS_HYDROMETEOR_entr_0D

  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (2D)
  subroutine ATMOS_HYDROMETEOR_entr_2D( &
       entr, &
       temp, &
       pres, &
       q,    &
       Rq    )
    use scale_const, only: &
       EPS   => CONST_EPS,   &
       PRE00 => CONST_PRE00, &
       TEM00 => CONST_TEM00, &
       CPdry => CONST_CPdry, &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       LHV0  => CONST_LHV0,  &
       LHF0  => CONST_LHF0,  &
       PSAT0 => CONST_PSAT0
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(out) :: entr(IA,JA)    !< entropy            [J/K]
    real(RP), intent(in)  :: temp(IA,JA)    !< temperature        [K]
    real(RP), intent(in)  :: pres(IA,JA)    !< pressure           [Pa]
    real(RP), intent(in)  :: q   (IA,JA,QA) !< mass concentration [kg/kg]
    real(RP), intent(in)  :: Rq  (QA)       !< gas constant       [J/kg/K]

    real(RP) :: qdry, qliq, qice, Rtot
    real(RP) :: logT_T0, pres_dry, pres_vap

    integer  :: i, j, iqw
    !---------------------------------------------------------------------------

    ! dry air + vapor
    do j = JSB, JEB
    do i = ISB, IEB

       logT_T0 = log( temp(i,j) / TEM00 )

       qdry = 1.0_RP
       Rtot = 0.0_RP
       do iqw = 1, QA
          qdry = qdry - q(i,j,iqw)
          Rtot = Rtot + q(i,j,iqw) * Rq(iqw)
       enddo
       Rtot = Rtot + Rdry * qdry

       ! dry air + vapor
       pres_dry = max( pres(i,j) * qdry    * Rdry / Rtot, EPS )
       entr(i,j) = qdry    * CPdry * logT_T0 &
                 - qdry    * Rdry  * log( pres_dry / PRE00 )

       if ( I_QV > 0 ) then
          pres_vap = max( pres(i,j) * q(i,j,I_QV) * Rvap / Rtot, EPS )
          entr(i,j) = entr(i,j) + q(i,j,I_QV) * CP_VAPOR * logT_T0 &
                                - q(i,j,I_QV) * Rvap  * log( pres_vap / PSAT0 ) &
                                + q(i,j,I_QV) * LHV0 / TEM00
       endif

       ! liquid water
       if ( QLS > 0 ) then
          do iqw = QLS, QLE
             entr(i,j) = entr(i,j) + q(i,j,iqw) * CP_WATER * logT_T0
          enddo
       endif

       ! ice
       if ( QIS > 0 ) then
          do iqw = QIS, QIE
             entr(i,j) = entr(i,j) + q(i,j,iqw) * CP_ICE * logT_T0 &
                                   - q(i,j,iqw) * LHF0 / TEM00
          enddo
       endif

    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_entr_2D

  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (3D)
  subroutine ATMOS_HYDROMETEOR_entr_3D( &
       entr, &
       temp, &
       pres, &
       q,    &
       Rq    )
    use scale_const, only: &
       EPS   => CONST_EPS,   &
       PRE00 => CONST_PRE00, &
       TEM00 => CONST_TEM00, &
       CPdry => CONST_CPdry, &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       LHV0  => CONST_LHV0,  &
       LHF0  => CONST_LHF0,  &
       PSAT0 => CONST_PSAT0
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(out) :: entr(KA,IA,JA)    !< entropy            [J/K]
    real(RP), intent(in)  :: temp(KA,IA,JA)    !< temperature        [K]
    real(RP), intent(in)  :: pres(KA,IA,JA)    !< pressure           [Pa]
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) !< mass concentration [kg/kg]
    real(RP), intent(in)  :: Rq  (QA)          !< gas constant       [J/kg/K]

    real(RP) :: qdry, qliq, qice, Rtot
    real(RP) :: logT_T0, pres_dry, pres_vap

    integer  :: k, i, j, iqw
    !---------------------------------------------------------------------------

    ! dry air + vapor
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE

       logT_T0 = log( temp(k,i,j) / TEM00 )

       qdry = 1.0_RP
       Rtot = 0.0_RP
       do iqw = 1, QA
          qdry = qdry - q(k,i,j,iqw)
          Rtot = Rtot + q(k,i,j,iqw) * Rq(iqw)
       enddo
       Rtot = Rtot + Rdry * qdry

       ! dry air + vapor
       pres_dry = max( pres(k,i,j) * qdry    * Rdry / Rtot, EPS )
       entr(k,i,j) = qdry    * CPdry * logT_T0 &
                   - qdry    * Rdry  * log( pres_dry / PRE00 )

       if ( I_QV > 0 ) then
          pres_vap = max( pres(k,i,j) * q(k,i,j,I_QV) * Rvap / Rtot, EPS )
          entr(k,i,j) = entr(k,i,j) + q(k,i,j,I_QV) * CP_VAPOR * logT_T0 &
                                    - q(k,i,j,I_QV) * Rvap  * log( pres_vap / PSAT0 ) &
                                    + q(k,i,j,I_QV) * LHV0 / TEM00
       endif

       ! liquid water
       if ( QLS > 0 ) then
          do iqw = QLS, QLE
             entr(k,i,j) = entr(k,i,j) + q(k,i,j,iqw) * CP_WATER * logT_T0
          enddo
       endif

       ! ice
       if ( QIS > 0 ) then
          do iqw = QIS, QIE
             entr(k,i,j) = entr(k,i,j) + qice * CP_ICE * logT_T0 &
                                       - qice * LHF0 / TEM00
          enddo
       endif

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROMETEOR_entr_3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_HYDROMETEOR_diagnose_number_concentration( &
       QTRC )
    use scale_const, only: &
       PI => CONST_PI
    implicit none

    real(RP), intent(inout) :: QTRC(:,:,:,:)

    real(RP), parameter :: Dc   =  20.D-6  ! typical particle diameter for cloud  [m]
    real(RP), parameter :: Dr   = 200.D-6  ! typical particle diameter for rain   [m]
    real(RP), parameter :: Di   =  80.D-6  ! typical particle diameter for ice    [m]
    real(RP), parameter :: Ds   =  80.D-6  ! typical particle diameter for snow   [m]
    real(RP), parameter :: Dg   = 200.D-6  ! typical particle diameter for grapel [m]
    real(RP), parameter :: RHOw = 1000.D0  ! typical density for warm particles   [kg/m3]
    real(RP), parameter :: RHOf =  100.D0  ! typical density for frozen particles [kg/m3]
    real(RP), parameter :: RHOg =  400.D0  ! typical density for grapel particles [kg/m3]
    real(RP), parameter :: b    =  3.D0    ! assume spherical form

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
