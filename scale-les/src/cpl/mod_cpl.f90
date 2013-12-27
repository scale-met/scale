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
  use mod_precision
  use mod_index
  use mod_tracer
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
       I_EMIT,             &
       I_ALB,              &
       I_TCS,              &
       I_DZg,              &
       I_Z0M,              &
       I_Z0H,              &
       I_Z0E,              &
       LAND_PROPERTY
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd, &
       CPL_vars_setup,             &
       CPL_vars_restart_read,      &
       CPL_vars_merge,             &
       CPL_vars_fillhalo,          &
       CPL_putAtm,                 &
       CPL_putLnd,                 &
       CPL_flushAtm,               &
       CPL_flushLnd,               &
       CPL_AtmLnd_flushCPL
    use mod_cpl_atmos_land, only: &
       CPL_AtmLnd_setup, &
       CPL_AtmLnd_unsolve
    implicit none
    !---------------------------------------------------------------------------

    call CPL_vars_setup
    call CPL_vars_restart_read

    call CPL_AtmLnd_setup

    if( sw_AtmLnd ) then
       call CPL_putATM( &
          DENS, MOMX, MOMY, MOMZ,    &
          RHOT, QTRC, PREC, SWD, LWD )
       call CPL_putLnd( TG           (:,:),        & ! [IN]
                        QvEfc        (:,:),        & ! [IN]
                        LAND_PROPERTY(:,:,I_EMIT), & ! [IN]
                        LAND_PROPERTY(:,:,I_ALB),  & ! [IN]
                        LAND_PROPERTY(:,:,I_TCS),  & ! [IN]
                        LAND_PROPERTY(:,:,I_DZg),  & ! [IN]
                        LAND_PROPERTY(:,:,I_Z0M),  & ! [IN]
                        LAND_PROPERTY(:,:,I_Z0H),  & ! [IN]
                        LAND_PROPERTY(:,:,I_Z0E)   ) ! [IN]
       call CPL_AtmLnd_unsolve
    endif

    call CPL_vars_merge

    call CPL_flushAtm
    call CPL_flushLnd
    call CPL_AtmLnd_flushCPL

    call CPL_vars_fillhalo

    return
  end subroutine CPL_setup

  !-----------------------------------------------------------------------------
  !> CPL calcuration
  subroutine CPL_calc
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd, &
       CPL_TYPE_AtmLnd,     &
       CPL_vars_fillhalo,   &
       CPL_vars_merge,      &
       CPL_vars_history,    &
       CPL_flushAtm,        &
       CPL_flushLnd,        &
       CPL_AtmLnd_flushCPL
    use mod_cpl_atmos_land, only: &
       CPL_AtmLnd_solve,   &
       CPL_AtmLnd_unsolve
    implicit none
    !---------------------------------------------------------------------------

    !########## Coupler Atoms-Land ##########
    call TIME_rapstart('CPL Atmos-Land')
    if( sw_AtmLnd ) then
       call CPL_AtmLnd_solve
    endif
    call TIME_rapend  ('CPL Atmos-Land')

    !########## merge Land-Ocean ##########
    call CPL_vars_merge

    !########## flush Coupler values ##########
    call CPL_flushAtm
    call CPL_flushLnd
    call CPL_AtmLnd_flushCPL

    !########## Fill HALO ##########
    call CPL_vars_fillhalo

    !########## History & Monitor ##########
    call TIME_rapstart('CPL History')
       call CPL_vars_history
    call TIME_rapend  ('CPL History')

    return
  end subroutine CPL_calc

end module mod_cpl
