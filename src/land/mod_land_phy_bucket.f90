!-------------------------------------------------------------------------------
!> module LAND / Physics Bucket
!!
!! @par Description
!!          bucket-type land physics module
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_land_phy_bucket
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_land.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_setup
  public :: LAND_PHY

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  real(RP), public, save :: ROFF   (IA,JA) ! run-off water [kg/m2]
  real(RP), public, save :: STRG   (IA,JA) ! water storage [kg/m2]
  real(RP), public, save :: STRGMAX(IA,JA) ! maximum water storage [kg/m2]
  real(RP), public, save :: STRGCRT(IA,JA) ! critical water storage [kg/m2]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

  ! limiter
  real(RP), private, parameter :: BETA_MIN = 0.0E-8_RP
  real(RP), private, parameter :: BETA_MAX = 1.0_RP

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_PHY_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_land_vars, only: &
       LAND_RESTART_IN_BASENAME, &
       LAND_TYPE_PHY
    implicit none

    logical  :: dummy

    NAMELIST / PARAM_LAND_BUCKET / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[BUCKET]/Categ[LAND]'

    if ( LAND_TYPE_PHY /= 'BUCKET' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx LAND_TYPE_PHY is not BUCKET. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_BUCKET,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_BUCKET. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LAND_BUCKET)

    ROFF   (:,:) = 0.0_RP
    STRG   (:,:) = 150.0_RP
    STRGMAX(:,:) = 150.0_RP
    STRGCRT(:,:) = STRGMAX * 0.75_RP

    return
  end subroutine LAND_PHY_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY
    use mod_const, only: &
      DWATR => CONST_DWATR, &
      CL    => CONST_CL
    use mod_time, only: &
      dt => TIME_DTSEC_LAND
    use mod_land_vars, only: &
      SFLX_GH,    &
      SFLX_PREC,  &
      SFLX_QVLnd, &
      TG,         &
      QvEfc,      &
      HCS,        &
      DZg
 
    implicit none

    integer :: i,j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE

      ! update water storage
      STRG(i,j)  = STRG(i,j) + ( SFLX_PREC(i,j) - SFLX_QVLnd(i,j) ) * dt

      if( STRG(i,j) > STRGMAX(i,j) ) then
        ROFF(i,j) = ROFF(i,j) + STRG(i,j) - STRGMAX(i,j)
        STRG(i,j) = STRGMAX(i,j)
      endif

      ! update moisture efficiency
      QvEfc(i,j) = BETA_MAX
      if( STRG(i,j) < STRGCRT(i,j) ) then
        QvEfc(i,j) = max( STRG(i,j)/STRGCRT(i,j), BETA_MIN )
      endif

      ! update ground temperature
      TG(i,j) = TG(i,j) - 2.0_RP * SFLX_GH(i,j) / ( ( 1.0_RP - STRGMAX(i,j) * 1.0E-3_RP ) * HCS(i,j) + STRG(i,j) * 1.0E-3_RP * DWATR * CL * DZg(i,j) ) * dt

    end do
    end do

    return
  end subroutine LAND_PHY

end module mod_land_phy_bucket
