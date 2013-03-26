!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          DYCOMS-II Parametarized Radiative heating
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-26 (H.Yashiro)  [new]
!! @li      2012-06-10 (Y.Miyamoto) [mod] bug-fix and some modification
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_setup
  public :: ATMOS_PHY_RD

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'
  include 'inc_precision.h'
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
  subroutine ATMOS_PHY_RD_setup
    return
  end subroutine ATMOS_PHY_RD_setup

  !-----------------------------------------------------------------------------
  ! Parametarized Radiative heating
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD
    use mod_time, only: &
       dtrd => TIME_DTSEC_ATMOS_PHY_RD
    use mod_grid, only : &
       CZ   => GRID_CZ
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       ATMOS_vars_total,   &
       DENS, &
       RHOT, &
       QTRC
    implicit none

    real(RP) :: TEMP_t  (KA,IA,JA) ! tendency rho*theta     [K*kg/m3/s]

    integer :: k, i, j

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Parametarized Radiation of RICO'

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
            RHOT(k,i,j) = RHOT(k,i,j) &
                        - dtrd * 2.5_RP/86400.0_RP * DENS(k,i,j)
            TEMP_t(k,i,j) = - 2.5_RP/86400.0_RP
       enddo

    enddo
    enddo

    ! fill KHALO
    do j  = JS, JE
    do i  = IS, IE
       RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
       RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
    enddo
    enddo
    ! fill IHALO & JHALO
    call COMM_vars8( RHOT(:,:,:), 5 )
    call COMM_wait ( RHOT(:,:,:), 5 )

    call ATMOS_vars_total

    call HIST_in( TEMP_t  (:,:,:), 'TEMP_t_rd', 'tendency of temp in rd', 'K*kg/m3/s', dtrd )

    return
  end subroutine ATMOS_PHY_RD

end module mod_atmos_phy_rd
