!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics Subset
!!
!! @par Description
!!          Subsets for Microphysics
!!          Seiki scheme
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-06 (H.Yashiro) [new] add dummy
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_mpsub
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MPsub_setup
  public :: ATMOS_PHY_MPsub_MP2RD
  public :: ATMOS_PHY_MPsub_CloudFraction
  public :: ATMOS_PHY_MPsub_EffectiveRadius

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public, parameter :: MP_QA = 1

  integer,  public, parameter :: I_mp_QC = 1

  integer,  public, save :: I_MP2ALL(MP_QA)
  data I_MP2ALL / -999 /

  real(RP), public, save :: MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]

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
  !> setup
  subroutine ATMOS_PHY_MPsub_setup
    use mod_const, only: &
       CONST_DWATR
    implicit none
    !---------------------------------------------------------------------------

    MP_DENS(:) = CONST_DWATR ! DUMMY!!!

    return
  end subroutine ATMOS_PHY_MPsub_setup

  !-----------------------------------------------------------------------------
  !> Make look-up table between hydrometeor tracer and particle type in radiation scheme
  subroutine ATMOS_PHY_MPsub_MP2RD( &
       I_MP2RD )
    implicit none

    integer, intent(out) :: I_MP2RD(MP_QA)
    !---------------------------------------------------------------------------

    I_MP2RD(:) = 1 ! DUMMY!!!

    return
  end subroutine ATMOS_PHY_MPsub_MP2RD

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MPsub_CloudFraction( &
       cldfrac, &
       QTRC     )
    use mod_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QA)

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.D0
       do iq = 1, MP_QA
          qhydro = qhydro + QTRC(k,i,j,I_MP2ALL(iq))
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-EPS)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MPsub_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MPsub_EffectiveRadius( &
       Re,  &
       QTRC )
    implicit none

    real(RP), intent(out) :: Re  (KA,IA,JA,MP_QA) ! effective radius
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    !---------------------------------------------------------------------------

    Re(:,:,:,:) = 8.E-6_RP ! DUMMY!!!

    return
  end subroutine ATMOS_PHY_MPsub_EffectiveRadius

end module mod_atmos_phy_mpsub 
