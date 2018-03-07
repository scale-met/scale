!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics / Cumulus Parameterization / Common
!!
!! @par Description
!!          Common module for Cumulus convection parameterization
!!          Running mean of vertical wind velocity
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_cp_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CP_common_setup
  public :: ATMOS_PHY_CP_common_wmean

!  interface ATMOS_PHY_CP_common_wmean
!     module procedure ATMOS_PHY_CP_common_wmean
!  end interface ATMOS_PHY_CP_common_wmean

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
  ! tuning parameter
  logical,  private :: PARAM_ATMOS_PHY_CP_wadapt = .true.
  integer,  private :: PARAM_ATMOS_PHY_CP_w_time = 16
  !------------------------------------------------------------------------------
contains
  !------------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CP_common_setup ()
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_CP_COMMON / &
       PARAM_ATMOS_PHY_CP_wadapt, &
       PARAM_ATMOS_PHY_CP_w_time

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CUMULUS] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** CP-COMMON'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_COMMON,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_CP_COMMON. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_CP_COMMON)

    ! output parameter lists
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) "*** Use running mean of w in adaptive timestep?     : ", PARAM_ATMOS_PHY_CP_wadapt
    if( IO_L ) write(IO_FID_LOG,*) "*** Fixed time scale for running mean of w          : ", PARAM_ATMOS_PHY_CP_w_time

    return
  end subroutine ATMOS_PHY_CP_common_setup

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_CP_wmean
  !! running mean vertical wind velocity
  !! comment for W0 imported from WRF
  !!
  !! TST IS THE NUMBER OF TIME STEPS IN 10 MINUTES...W0AVG IS CLOSE TO A
  !! RUNNING MEAN VERTICAL VELOCITY...NOTE THAT IF YOU CHANGE TST, IT WIL
  !! CHANGE THE FREQUENCY OF THE CONVECTIVE INTITIATION CHECK (SEE BELOW)
  !! NOTE THAT THE ORDERING OF VERTICAL LAYERS MUST BE REVERSED FOR W0AVG
  !< BECAUSE THE ORDERING IS REVERSED IN KFPARA.
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_CP_common_wmean( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS,    &
       MOMZ,    &
       W0_mean  )
    use scale_time , only :&
       TIME_DTSEC,             &
       CP_DTSEC => TIME_DTSEC_ATMOS_PHY_CP
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)    :: DENS   (KA,IA,JA)
    real(RP), intent(in)    :: MOMZ   (KA,IA,JA)
    real(RP), intent(inout) :: W0_mean(KA,IA,JA)

    real(RP) :: W0
    real(RP) :: fact_W0_mean, fact_W0

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( PARAM_ATMOS_PHY_CP_wadapt ) then
       fact_W0_mean = 2.0_RP * max(CP_DTSEC,TIME_DTSEC) - TIME_DTSEC
       fact_W0      = TIME_DTSEC
    else ! w_time is tuning parameter
       fact_W0_mean = real(PARAM_ATMOS_PHY_CP_w_time,RP)
       fact_W0      = 1.0_RP
    endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       W0 = 0.5_RP * ( MOMZ(k,i,j) + MOMZ(k-1,i,j) ) / DENS(k,i,j)

       W0_mean(k,i,j) = ( W0_mean(k,i,j) * fact_W0_mean &
                        + W0             * fact_W0      ) / ( fact_W0_mean + fact_W0 )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_CP_common_wmean

end module scale_atmos_phy_cp_common
