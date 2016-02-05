!-----------------------------------------------------------------------------
!> Runge-Kutta loop
!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          Runge-Kutta for Atmospheric dynamical process
!!          NONE
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-09-30 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_rk_none
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
  public :: ATMOS_DYN_rk_none_regist
  public :: ATMOS_DYN_rk_none_setup
  public :: ATMOS_DYN_rk_none

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
  integer, parameter :: VA_NONE = 0

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_rk_none_regist( &
       ATMOS_DYN_TYPE, &
       VA_out, &
       VAR_NAME, VAR_DESC, VAR_UNIT )
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    character(len=*),       intent(in)  :: ATMOS_DYN_TYPE
    integer,                intent(out) :: VA_out !< number of prognostic variables
    character(len=H_SHORT), intent(out) :: VAR_NAME(:) !< name  of the variables
    character(len=H_MID),   intent(out) :: VAR_DESC(:) !< desc. of the variables
    character(len=H_SHORT), intent(out) :: VAR_UNIT(:) !< unit  of the variables

    if( IO_L ) write(IO_FID_LOG,*) '*** Dyn NONE Register'

    if ( ATMOS_DYN_TYPE .ne. 'NONE' ) then
       write(*,*) 'xxx ATMOS_DYN_TYPE is not NONE. Check!'
       call PRC_MPIstop
    endif

    VA_out = VA_NONE

    return
  end subroutine ATMOS_DYN_rk_none_regist

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_rk_none_setup( &
       BND_W, &
       BND_E, &
       BND_S, &
       BND_N  )
    implicit none

    logical, intent(in) :: BND_W
    logical, intent(in) :: BND_E
    logical, intent(in) :: BND_S
    logical, intent(in) :: BND_N
    !---------------------------------------------------------------------------

    return
  end subroutine ATMOS_DYN_rk_none_setup

  !-----------------------------------------------------------------------------
  !> Runge-Kutta loop
  subroutine ATMOS_DYN_rk_none( &
       DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
       PROG_RK,                                     &
       mflx_hi, tflx_hi,                            &
       DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
       DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
       PROG0, PROG,                                 &
       Rtot, CVtot, CORIOLI,                        &
       num_diff, divdmp_coef,                       &
       FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T, &
       FLAG_FCT_ALONG_STREAM,                       &
       CDZ, FDZ, FDX, FDY,                          &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
       PHI, GSQRT, J13G, J23G, J33G, MAPF,          &
       REF_pres, REF_dens,                          &
       BND_W, BND_E, BND_S, BND_N,                  &
       dtrk, dt                                     )
    use scale_grid_index
    use scale_index
    implicit none

    real(RP), intent(out) :: DENS_RK(KA,IA,JA)   ! prognostic variables
    real(RP), intent(out) :: MOMZ_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMX_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMY_RK(KA,IA,JA)   !
    real(RP), intent(out) :: RHOT_RK(KA,IA,JA)   !

    real(RP), intent(out) :: PROG_RK(KA,IA,JA,VA)  !

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

    logical,  intent(in)  :: FLAG_FCT_RHO
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

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    ! DUMMY

    return
  end subroutine ATMOS_DYN_rk_none

end module scale_atmos_dyn_rk_none
