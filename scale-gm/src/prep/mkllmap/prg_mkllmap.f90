!-------------------------------------------------------------------------------
!
!+  Program mkllmap
!
!-------------------------------------------------------------------------------
program prg_mkllmap
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !
  !++ Current Corresponding Author :
  !
  !++ History:
  !       This program originate from 'cnvlatlon.f90'(by H.Tomita) ver 4.38 .
  !
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !                 11-11-09  H.Yashiro [mod] Avoid arc-cos, precise calculation
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules (shared)
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID,     &
     ADM_MULTI_PRC,   &
     ADM_proc_init,   &
     ADM_proc_stop,   &
     ADM_proc_finish, &
     ADM_setup,       &
     ADM_CTL_FID,     &
     ADM_MAXFNAME
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_cnst, only: &
     CNST_setup
  use mod_grd, only: &
     GRD_setup
  use mod_latlon, only: &
     LATLON_setup, &
     LATLON_ico_setup
  implicit none

  character(len=ADM_MAXFNAME) :: output_dir   = './'
  logical                     :: use_quadprec = .false.

  namelist /mkllmap_param/ &
      use_quadprec, &
      output_dir

  integer :: ierr
  !=============================================================================
  !
  !--- start process
  !
  call ADM_proc_init(ADM_MULTI_PRC)
  !
  !--- < admin module setup > ---
  call ADM_setup('mkllmap.cnf')
  !
  call FIO_setup
  !
  !--- < comm module setup > ---
  call COMM_setup
  !
  !--- < cnst module setup > ---
  call CNST_setup
  !
  !--- < grid module setup > ---
  call GRD_setup

  !--- read parameters
  write(ADM_LOG_FID,*)
  write(ADM_LOG_FID,*) '+++ Program[mkllmap]/Category[tool]'
  rewind(ADM_CTL_FID)
  read(ADM_CTL_FID,nml=MKLLMAP_PARAM,iostat=ierr)
  if ( ierr < 0 ) then
     write(ADM_LOG_FID,*) '*** MKLLMAP_PARAM is not specified. use default.'
  elseif( ierr > 0 ) then
     write(*,          *) 'xxx Not appropriate names in namelist MKLLMAP_PARAM. STOP.'
     write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist MKLLMAP_PARAM. STOP.'
     call ADM_proc_stop
  endif
  write(ADM_LOG_FID,nml=MKLLMAP_PARAM)

  call LATLON_ico_setup

  if ( use_quadprec ) then
     call LATLON_setup( output_dir, 'mkllmap_q' )
  else
     call LATLON_setup( output_dir, 'mkllmap' )
  endif


  call ADM_proc_finish

end program prg_mkllmap
!-------------------------------------------------------------------------------
