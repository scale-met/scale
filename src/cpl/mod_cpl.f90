!-------------------------------------------------------------------------------
!> module CPL driver
!!
!! @par Description
!!          CPL module driver
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_cpl
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
  !++ Public procedure
  !
  public :: CPL_setup
  public :: CPL_calc

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  !> Setup
  subroutine CPL_setup
    use mod_atmos_vars, only: &
       DENS, & 
       MOMX, & 
       MOMY, & 
       MOMZ, & 
       RHOT, & 
       QTRC
    use mod_atmos_vars_sf, only: &
       PREC, &
       SWD,  &
       LWD
    use mod_land_vars, only: &
       TG,    &
       QvEfc, &
       EMIT,  &
       ALB, TCS, DZg, &
       Z0M, Z0H, Z0E
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd, &
       CPL_vars_setup,             &
       CPL_vars_restart_read,      &
       CPL_vars_fillhalo
    use mod_cpl_atmos_land, only: &
       CPL_AtmLnd_setup,   &
       CPL_AtmLnd_unsolve, &
       CPL_AtmLnd_putAtm,  &
       CPL_AtmLnd_putLnd
    implicit none
    !---------------------------------------------------------------------------

    call CPL_vars_setup

    call CPL_vars_restart_read

    call CPL_AtmLnd_setup

    if( sw_AtmLnd ) then
       call CPL_AtmLnd_putATM( &
          DENS, MOMX, MOMY, MOMZ,    &
          RHOT, QTRC, PREC, SWD, LWD )
       call CPL_AtmLnd_putLnd( &
          TG, QvEfc, EMIT, &
          ALB, TCS, DZg,   &
          Z0M, Z0H, Z0E    )
       call CPL_AtmLnd_unsolve
    endif

    call CPL_vars_fillhalo

    return
  end subroutine CPL_setup

  !-----------------------------------------------------------------------------
  !> CPL calcuration
  subroutine CPL_calc
    use mod_time, only: &
       do_AtmLnd => TIME_DOCPL_AtmLnd
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd, &
       CPL_TYPE_AtmLnd, &
       CPL_vars_fillhalo
    use mod_cpl_atmos_land, only: &
       CPL_AtmLnd_solve,   &
       CPL_AtmLnd_unsolve
    implicit none
    !---------------------------------------------------------------------------

    !########## Coupler Atoms-Land ##########
    call TIME_rapstart('CPL Atmos-Land')
    if( sw_AtmLnd ) then
       if( do_AtmLnd ) call CPL_AtmLnd_solve
    endif
    call TIME_rapend  ('CPL Atmos-Land')

    !########## Fill HALO ##########
    call CPL_vars_fillhalo

    return
  end subroutine CPL_calc

end module mod_cpl
