!-------------------------------------------------------------------------------
!> module ocean / physics / surface albedo
!!
!! @par Description
!!          Ocean surface albedo common module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_phy_albedo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_PHY_ALBEDO_const_setup
  public :: OCEAN_PHY_ALBEDO_seaice_setup
  public :: OCEAN_PHY_ALBEDO_const
  public :: OCEAN_PHY_ALBEDO_seaice

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
  real(RP), private :: OCEAN_PHY_ALBEDO_IR_dir         = 0.05_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_IR_dif         = 0.05_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_NIR_dir        = 0.07_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_NIR_dif        = 0.06_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_VIS_dir        = 0.07_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_VIS_dif        = 0.06_RP

  real(RP), private :: OCEAN_PHY_ALBEDO_seaice_IR_dir  = 0.05_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_seaice_IR_dif  = 0.05_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_seaice_NIR_dir = 0.60_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_seaice_NIR_dif = 0.60_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_seaice_VIS_dir = 0.80_RP
  real(RP), private :: OCEAN_PHY_ALBEDO_seaice_VIS_dif = 0.80_RP

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ALBEDO_const_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_OCEAN_PHY_ALBEDO_const / &
       OCEAN_PHY_ALBEDO_IR_dir,  &
       OCEAN_PHY_ALBEDO_IR_dif,  &
       OCEAN_PHY_ALBEDO_NIR_dir, &
       OCEAN_PHY_ALBEDO_NIR_dif, &
       OCEAN_PHY_ALBEDO_VIS_dir, &
       OCEAN_PHY_ALBEDO_VIS_dif

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ALBEDO_const_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ALBEDO_const,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ALBEDO_const_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ALBEDO_const_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ALBEDO_const. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ALBEDO_const)

    return
  end subroutine OCEAN_PHY_ALBEDO_const_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ALBEDO_seaice_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_OCEAN_PHY_ALBEDO_seaice / &
       OCEAN_PHY_ALBEDO_seaice_IR_dir,  &
       OCEAN_PHY_ALBEDO_seaice_IR_dif,  &
       OCEAN_PHY_ALBEDO_seaice_NIR_dir, &
       OCEAN_PHY_ALBEDO_seaice_NIR_dif, &
       OCEAN_PHY_ALBEDO_seaice_VIS_dir, &
       OCEAN_PHY_ALBEDO_seaice_VIS_dif

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ALBEDO_seaice_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ALBEDO_seaice,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ALBEDO_seaice_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ALBEDO_seaice_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ALBEDO_seaice. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ALBEDO_seaice)

    return
  end subroutine OCEAN_PHY_ALBEDO_seaice_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ALBEDO_const( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       SFC_albedo     )
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(out) :: SFC_albedo(OIA,OJA,N_RAD_DIR,N_RAD_RGN) ! surface albedo (0-1)

    integer :: i, j
    !---------------------------------------------------------------------------
    !$acc data copyout(SFC_albedo)

    !$acc kernels
    do j = OJS, OJE
    do i = OIS, OIE
       SFC_albedo(i,j,I_R_direct ,I_R_IR ) = OCEAN_PHY_ALBEDO_IR_dir
       SFC_albedo(i,j,I_R_diffuse,I_R_IR ) = OCEAN_PHY_ALBEDO_IR_dif
       SFC_albedo(i,j,I_R_direct ,I_R_NIR) = OCEAN_PHY_ALBEDO_NIR_dir
       SFC_albedo(i,j,I_R_diffuse,I_R_NIR) = OCEAN_PHY_ALBEDO_NIR_dif
       SFC_albedo(i,j,I_R_direct ,I_R_VIS) = OCEAN_PHY_ALBEDO_VIS_dir
       SFC_albedo(i,j,I_R_diffuse,I_R_VIS) = OCEAN_PHY_ALBEDO_VIS_dif
    enddo
    enddo
    !$acc end kernels

    !$acc end data
    return
  end subroutine OCEAN_PHY_ALBEDO_const

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ALBEDO_seaice( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       SFC_albedo     )
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(out) :: SFC_albedo(OIA,OJA,N_RAD_DIR,N_RAD_RGN) ! surface albedo (0-1)

    integer :: i, j
    !---------------------------------------------------------------------------
    !$acc data copyout(SFC_albedo)

    !$acc kernels
    do j = OJS, OJE
    do i = OIS, OIE
       SFC_albedo(i,j,I_R_direct ,I_R_IR ) = OCEAN_PHY_ALBEDO_seaice_IR_dir
       SFC_albedo(i,j,I_R_diffuse,I_R_IR ) = OCEAN_PHY_ALBEDO_seaice_IR_dif
       SFC_albedo(i,j,I_R_direct ,I_R_NIR) = OCEAN_PHY_ALBEDO_seaice_NIR_dir
       SFC_albedo(i,j,I_R_diffuse,I_R_NIR) = OCEAN_PHY_ALBEDO_seaice_NIR_dif
       SFC_albedo(i,j,I_R_direct ,I_R_VIS) = OCEAN_PHY_ALBEDO_seaice_VIS_dir
       SFC_albedo(i,j,I_R_diffuse,I_R_VIS) = OCEAN_PHY_ALBEDO_seaice_VIS_dif
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    return
  end subroutine OCEAN_PHY_ALBEDO_seaice

end module scale_ocean_phy_albedo
