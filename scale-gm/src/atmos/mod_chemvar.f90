!-------------------------------------------------------------------------------
!> Module tracer variable
!!
!! @par Description
!!          This module contains the chemical or general-perpose tracer variables
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_chemvar
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CHEMVAR_setup
  public :: chemvar_getid

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,                 public, parameter   :: CHEM_TRC_vlim = 100
  integer,                 public              :: CHEM_TRC_vmax = 1
  character(len=H_SHORT),  public, allocatable :: CHEM_TRC_name(:) ! short name  of tracer
  character(len=H_MID),    public, allocatable :: CHEM_TRC_desc(:) ! description of tracer

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CHEMVAR_setup
    use mod_adm, only: &
       ADM_CTL_FID, &
       ADM_proc_stop
    implicit none

    namelist /CHEMVARPARAM/ &
       CHEM_TRC_vmax

    integer :: nq
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[chemvar]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=CHEMVARPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** CHEMVARPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist CHEMVARPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist CHEMVARPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=CHEMVARPARAM)

    allocate( CHEM_TRC_name(CHEM_TRC_vmax) )
    allocate( CHEM_TRC_desc(CHEM_TRC_vmax) )

    do nq = 1, CHEM_TRC_vmax
       write(CHEM_TRC_name(nq),'(A,I3.3)') 'passive', nq
       write(CHEM_TRC_desc(nq),'(A,I3.3)') 'passive_tracer_no', nq
    enddo

    return
  end subroutine CHEMVAR_setup

  !-----------------------------------------------------------------------------
  function chemvar_getid( tracername )
    use mod_adm, only: &
       ADM_proc_stop
    implicit none

    character(len=*), intent(in) :: tracername
    integer                      :: chemvar_getid

    character(len=H_SHORT) :: tname
    integer           :: itrc
    !---------------------------------------------------------------------------

    tname = trim(tracername)

    chemvar_getid = -1

    do itrc = 1, CHEM_TRC_vmax
       if ( tname == CHEM_TRC_name(itrc) ) then
          chemvar_getid = itrc
          return
       endif
    enddo

    if ( chemvar_getid <= 0 ) then
       write(ADM_LOG_FID,*) 'xxx INDEX does not exist =>', tname
       call ADM_proc_stop
    endif

  end function chemvar_getid

end module mod_chemvar
