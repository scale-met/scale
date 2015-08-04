!-------------------------------------------------------------------------------
!> Program mkhgrid
!!
!! @par Description
!!          Making horizontal grid systems based on the icosahedral grid configuration
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program mkhgrid
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_LOG_FID,     &
     ADM_MULTI_PRC,   &
     ADM_proc_init,   &
     ADM_proc_finish, &
     ADM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use scale_const, only: &
     CONST_setup, &
     RADIUS => CONST_RADIUS
  use mod_grd, only: &
     GRD_input_hgrid,  &
     GRD_output_hgrid, &
     GRD_scaling
  use mod_gmtr, only: &
     GMTR_setup
  use mod_mkgrd, only: &
     MKGRD_setup,        &
     MKGRD_prerotate,    &
     MKGRD_stretch,      &
     MKGRD_shrink,       &
     MKGRD_rotate,       &
     MKGRD_gravcenter,   &
     MKGRD_diagnosis,    &
     MKGRD_IN_BASENAME,  &
     MKGRD_IN_io_mode,   &
     MKGRD_OUT_BASENAME, &
     MKGRD_OUT_io_mode
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  call ADM_proc_init(ADM_MULTI_PRC)

  !---< admin module setup >---
  call ADM_setup('mkhgrid.cnf')

  !---< I/O module setup >---
  call FIO_setup

  !---< comm module setup >---
  call COMM_setup

  !---< cnst module setup >---
  call CONST_setup

  !---< mkgrid module setup >---
  call MKGRD_setup

  !########## main ##########

  call GRD_input_hgrid( basename     = MKGRD_IN_BASENAME, &
                        input_vertex = .false.,           &
                        io_mode      = MKGRD_IN_io_mode   )

  call MKGRD_prerotate

  call MKGRD_stretch

  call MKGRD_shrink

  call MKGRD_rotate

  call MKGRD_gravcenter

  call GRD_output_hgrid( basename      = MKGRD_OUT_BASENAME, &
                         output_vertex = .true.,             &
                         io_mode       = MKGRD_OUT_io_mode   )

  !---< gmtr module setup >---
  call GRD_scaling( RADIUS )

  call GMTR_setup

  call MKGRD_diagnosis

  !########## Finalize ##########

  !--- all processes stop
  call ADM_proc_finish

  stop
end program mkhgrid
!-------------------------------------------------------------------------------
