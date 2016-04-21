!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamical scheme
!!
!! @par Description
!!          Dynamical scheme selecter for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-04-18 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tstep
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tstep_regist

  abstract interface
     !> setup
     subroutine setup( &
          scheme )
       character(len=*), intent(in) :: scheme
     end subroutine setup

     !> calculation values at next temporal step
     subroutine tstep( DENS_new, MOMZ_new, MOMX_new, MOMY_new, RHOT_new, & ! (out)
                       PROG_new,                                         & ! (out)
                       mflx_hi,  tflx_hi,                                & ! (out)
                       DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                       DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (in)
                       DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                       PROG0, PROG,                                      & ! (in)
                       Rtot, CVtot, CORIOLI,                             & ! (in)
                       num_diff, divdmp_coef, DDIV,                      & ! (in)
                       FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! (in)
                       FLAG_FCT_ALONG_STREAM,                            & ! (in)
                       CDZ, FDZ, FDX, FDY,                               & ! (in)
                       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                       PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! (in)
                       REF_pres, REF_dens,                               & ! (in)
                       BND_W, BND_E, BND_S, BND_N,                       & ! (in)
                       dtrk, dt                                          ) ! (in)
       use scale_precision
       use scale_grid_index
       use scale_index
       real(RP), intent(out) :: DENS_new(KA,IA,JA)   ! prognostic variables
       real(RP), intent(out) :: MOMZ_new(KA,IA,JA)   !
       real(RP), intent(out) :: MOMX_new(KA,IA,JA)   !
       real(RP), intent(out) :: MOMY_new(KA,IA,JA)   !
       real(RP), intent(out) :: RHOT_new(KA,IA,JA)   !
       real(RP), intent(out) :: PROG_new(KA,IA,JA,VA)  !

       real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3) ! mass flux
       real(RP), intent(out)   :: tflx_hi(KA,IA,JA,3) ! internal energy flux

       real(RP), intent(in),target :: DENS0(KA,IA,JA) ! prognostic variables at previous dynamical time step
       real(RP), intent(in),target :: MOMZ0(KA,IA,JA) !
       real(RP), intent(in),target :: MOMX0(KA,IA,JA) !
       real(RP), intent(in),target :: MOMY0(KA,IA,JA) !
       real(RP), intent(in),target :: RHOT0(KA,IA,JA) !

       real(RP), intent(in)  :: DENS(KA,IA,JA)      ! prognostic variables at previous RK step
       real(RP), intent(in)  :: MOMZ(KA,IA,JA)      !
       real(RP), intent(in)  :: MOMX(KA,IA,JA)      !
       real(RP), intent(in)  :: MOMY(KA,IA,JA)      !
       real(RP), intent(in)  :: RHOT(KA,IA,JA)      !

       real(RP), intent(in)  :: DENS_t(KA,IA,JA)    ! tendency
       real(RP), intent(in)  :: MOMZ_t(KA,IA,JA)    !
       real(RP), intent(in)  :: MOMX_t(KA,IA,JA)    !
       real(RP), intent(in)  :: MOMY_t(KA,IA,JA)    !
       real(RP), intent(in)  :: RHOT_t(KA,IA,JA)    !

       real(RP), intent(in)  :: PROG0(KA,IA,JA,VA)
       real(RP), intent(in)  :: PROG (KA,IA,JA,VA)

       real(RP), intent(in)  :: Rtot    (KA,IA,JA)  ! total R
       real(RP), intent(in)  :: CVtot   (KA,IA,JA)  ! total CV
       real(RP), intent(in)  :: CORIOLI (1, IA,JA)
       real(RP), intent(in)  :: num_diff(KA,IA,JA,5,3)
       real(RP), intent(in)  :: divdmp_coef
       real(RP), intent(in)  :: DDIV(KA,IA,JA)

       logical,  intent(in)  :: FLAG_FCT_MOMENTUM
       logical,  intent(in)  :: FLAG_FCT_T
       logical,  intent(in)  :: FLAG_FCT_ALONG_STREAM

       real(RP), intent(in)  :: CDZ (KA)
       real(RP), intent(in)  :: FDZ (KA-1)
       real(RP), intent(in)  :: FDX (IA-1)
       real(RP), intent(in)  :: FDY (JA-1)
       real(RP), intent(in)  :: RCDZ(KA)
       real(RP), intent(in)  :: RCDX(IA)
       real(RP), intent(in)  :: RCDY(JA)
       real(RP), intent(in)  :: RFDZ(KA-1)
       real(RP), intent(in)  :: RFDX(IA-1)
       real(RP), intent(in)  :: RFDY(JA-1)

       real(RP), intent(in)  :: PHI     (KA,IA,JA)   !< geopotential
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
       real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
       real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
       real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
       real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
       real(RP), intent(in)  :: REF_pres(KA,IA,JA)   !< reference pressure
       real(RP), intent(in)  :: REF_dens(KA,IA,JA)   !< reference density

       logical,  intent(in)  :: BND_W
       logical,  intent(in)  :: BND_E
       logical,  intent(in)  :: BND_S
       logical,  intent(in)  :: BND_N

       real(RP), intent(in)  :: dtrk
       real(RP), intent(in)  :: dt
     end subroutine tstep
  end interface

  procedure(setup), pointer :: ATMOS_DYN_Tstep_setup => NULL()
  public :: ATMOS_DYN_Tstep_setup
  procedure(tstep), pointer :: ATMOS_DYN_Tstep => NULL()
  public :: ATMOS_DYN_Tstep

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_Tstep_regist( &
       ATMOS_DYN_TYPE, &
       VA_out, &
       VAR_NAME, VAR_DESC, VAR_UNIT )
    use scale_precision
    use scale_grid_index
    use scale_index
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef DYNAMICS
    use NAME(scale_atmos_dyn_tstep_, DYNAMICS,), only: &
       NAME(ATMOS_DYN_rk_tstep_, DYNAMICS, _regist), &
       NAME(ATMOS_DYN_rk_tstep_, DYNAMICS, _setup), &
       NAME(ATMOS_DYN_rk_tstep_, DYNAMICS,)
#else
    use scale_atmos_dyn_tstep_fvm_heve, only: &
       ATMOS_DYN_Tstep_fvm_heve_regist, &
       ATMOS_DYN_Tstep_fvm_heve_setup, &
       ATMOS_DYN_Tstep_fvm_heve
    use scale_atmos_dyn_tstep_fvm_hevi, only: &
       ATMOS_DYN_Tstep_fvm_hevi_regist, &
       ATMOS_DYN_Tstep_fvm_hevi_setup, &
       ATMOS_DYN_Tstep_fvm_hevi
    use scale_atmos_dyn_tstep_fvm_hivi, only: &
       ATMOS_DYN_Tstep_fvm_hivi_regist, &
       ATMOS_DYN_Tstep_fvm_hivi_setup, &
       ATMOS_DYN_Tstep_fvm_hivi
#endif
    implicit none
    character(len=*),       intent(in)  :: ATMOS_DYN_TYPE
    integer,                intent(out) :: VA_out !< number of prognostic variables
    character(len=H_SHORT), intent(out) :: VAR_NAME(:) !< name  of the variables
    character(len=H_MID),   intent(out) :: VAR_DESC(:) !< desc. of the variables
    character(len=H_SHORT), intent(out) :: VAR_UNIT(:) !< unit  of the variables
    !---------------------------------------------------------------------------

#ifdef DYNAMICS
    NAME(ATMOS_DYN_Tstep_, DYNAMICS, _regist)( &
            ATMOS_DYN_TYPE )
    ATMOS_DYN_Tstep => NAME(ATMOS_DYN_Tstep_, DYNAMICS,)
    ATMOS_DYN_Tstep_setup => NAME(ATMOS_DYN_Tstep_, DYNAMICS, _setup)
#else
    select case ( ATMOS_DYN_TYPE )
    case ( 'FVM-HEVE', 'HEVE' )
       call ATMOS_DYN_Tstep_fvm_heve_regist( &
            ATMOS_DYN_TYPE, &
            VA_out, &
            VAR_NAME, VAR_DESC, VAR_UNIT )
       ATMOS_DYN_Tstep_setup => ATMOS_DYN_Tstep_fvm_heve_setup
       ATMOS_DYN_Tstep => ATMOS_DYN_Tstep_fvm_heve
    case ( 'FVM-HEVI', 'HEVI' )
       call ATMOS_DYN_Tstep_fvm_hevi_regist( &
            ATMOS_DYN_TYPE, &
            VA_out, &
            VAR_NAME, VAR_DESC, VAR_UNIT )
       ATMOS_DYN_Tstep_setup => ATMOS_DYN_Tstep_fvm_hevi_setup
       ATMOS_DYN_Tstep => ATMOS_DYN_Tstep_fvm_hevi
    case ( 'FVM-HIVI', 'HIVI' )
       call ATMOS_DYN_Tstep_fvm_hivi_regist( &
            ATMOS_DYN_TYPE, &
            VA_out, &
            VAR_NAME, VAR_DESC, VAR_UNIT )
       ATMOS_DYN_Tstep_setup => ATMOS_DYN_Tstep_fvm_hivi_setup
       ATMOS_DYN_Tstep => ATMOS_DYN_Tstep_fvm_hivi
    case ( 'OFF', 'NONE' )
       ! do nothing
    case default
       write(*,*) 'xxx ATMOS_DYN_TYPE is invalid: ', ATMOS_DYN_TYPE
       call PRC_MPIstop
    end select
#endif

    return
  end subroutine ATMOS_DYN_Tstep_regist

end module scale_atmos_dyn_tstep
