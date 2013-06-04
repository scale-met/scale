!-------------------------------------------------------------------------------
!> module OPTION
!!
!! @par Description
!!          Read argument and convert to namelist format
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
module mod_option
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ public procedure
  !
  public :: OPT_convert

  !-----------------------------------------------------------------------------
  !
  !++ public param & variable
  !
  integer, public :: OPT_fid !< fileunit number for namelist

  !-----------------------------------------------------------------------------
  !
  !++ private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ private param & variable
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Read argument and convert to namelist format
  subroutine OPT_convert( ninfile )
    implicit none

    integer, intent(out) :: ninfile !< nuber of input files

    character(LEN=128) :: argstr
    character(LEN=2)   :: snf

    integer :: command_argument_count

    integer :: n, ls, eq
    !---------------------------------------------------------------------------

    ninfile = 0

    write(OPT_fid,'(A)') '&OPTION'

    do n = 1, command_argument_count()
       call get_command_argument(n,argstr)
       ls = len_trim(argstr)

       if ( argstr(1:1) == '-' ) then
          if ( argstr(2:2) == '-' ) then                               ! '--option' format
             write(OPT_fid,'(A)') ' '//argstr(3:ls)//'=F'
          elseif ( argstr(2:3) == 'no' .OR. argstr(2:3) == 'NO' ) then ! '-nooption'/'-NOOPTION' format
             write(OPT_fid,'(A)') ' '//argstr(4:ls)//'=F'
          else                                                         ! '-option' format
             write(OPT_fid,'(A)') ' '//argstr(2:ls)//'=T'
          endif
       elseif( index(argstr,'=') == 0 ) then                           ! no '=' is filename
          ninfile = ninfile + 1
          write(snf,'(I2.2)') ninfile
          write(OPT_fid,'(A)') ' infile('//snf//')="'//trim(adjustl(argstr))//'"'
       elseif(      index(argstr,'/') /= 0 &
               .OR. index(argstr,' ') /= 0 ) then                      ! append "string"
          eq = index(argstr,'=')
          write(OPT_fid,'(A)') ' '//argstr(1:eq)//'"'//argstr(eq+1:ls)//'"'
       else                                                            ! without any change
          write(OPT_fid,'(A)') ' '//argstr(1:ls)
       endif
    enddo

    write(OPT_fid,'(A)') ' &END'
    rewind(OPT_fid)

    return
  end subroutine OPT_convert

end module mod_option
!-------------------------------------------------------------------------------
