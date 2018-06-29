!-------------------------------------------------------------------------------
!> Module tracer variable
!!
!! @par Description
!!         This module contains the chemical or general-perpose tracer variables
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
  use scale_io
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
  integer, public :: NCHEM_MAX =  0
  integer, public :: NCHEM_STR = -1
  integer, public :: NCHEM_END = -1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer                             :: CHEM_TRC_vmax = 0
  character(len=H_SHORT), allocatable :: CHEM_TRC_name(:) ! short name  of tracer
  character(len=H_MID),   allocatable :: CHEM_TRC_desc(:) ! description of tracer
  character(len=H_SHORT), allocatable :: CHEM_TRC_unit(:)


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CHEMVAR_setup
    use scale_prc, only: &
       PRC_abort
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    namelist /CHEMVARPARAM/ &
       CHEM_TRC_vmax

    integer :: nq
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[chemvar]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=CHEMVARPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** CHEMVARPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist CHEMVARPARAM. STOP.'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=CHEMVARPARAM)

    allocate( CHEM_TRC_name(max(CHEM_TRC_vmax,1)) )
    allocate( CHEM_TRC_desc(max(CHEM_TRC_vmax,1)) )
    allocate( CHEM_TRC_unit(max(CHEM_TRC_vmax,1)) )

    do nq = 1, CHEM_TRC_vmax
       write(CHEM_TRC_name(nq),'(A,I3.3)') 'passive', nq
       write(CHEM_TRC_desc(nq),'(A,I3.3)') 'passive_tracer_no', nq
       CHEM_TRC_unit(nq) = 'kg/kg'
    enddo

    if ( CHEM_TRC_vmax > 0 ) then
       call TRACER_regist( NCHEM_STR,                                           & ! [OUT]
                           CHEM_TRC_vmax,                                       & ! [IN]
                           CHEM_TRC_name(:), CHEM_TRC_desc(:), CHEM_TRC_unit(:) ) ! [IN]
       NCHEM_END = NCHEM_STR + CHEM_TRC_vmax - 1
       NCHEM_MAX = CHEM_TRC_vmax
    end if

    return
  end subroutine CHEMVAR_setup

  !-----------------------------------------------------------------------------
  function chemvar_getid( tracername )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: tracername
    integer                      :: chemvar_getid

    character(len=H_SHORT) :: tname
    integer                :: itrc
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
       write(*,*) 'xxx [chemvar_getid] INDEX does not exist =>', tname
       call PRC_abort
    endif

  end function chemvar_getid

end module mod_chemvar
