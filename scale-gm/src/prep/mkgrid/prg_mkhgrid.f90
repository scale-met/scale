!-------------------------------------------------------------------------------
!>
!! Program mkhgrid
!!
!! @par Description
!!          Making grid systems based on the icosahedral grid configuration
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2013-05-1  (H.Yashiro) NICAM-DC
!!
!<
program mkhgrid
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
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
  use mod_cnst, only: &
     CNST_setup, &
     RADIUS => CNST_ERADIUS
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
  call CNST_setup

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
