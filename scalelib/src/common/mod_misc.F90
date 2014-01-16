!-------------------------------------------------------------------------------
!> module MISC
!!
!! @par Description
!!          Miscellaneous tool module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-22 (H.Yashiro) [new]
!!
!<
module mod_misc
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MISC_valcheck

  interface MISC_valcheck
     module procedure MISC_valcheck_1D
     module procedure MISC_valcheck_2D
     module procedure MISC_valcheck_3D
  end interface MISC_valcheck

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
  !> Nan & extreme value checker (1D)
  subroutine MISC_valcheck_1D( &
       var,    &
       valmin, &
       valmax, &
       varname )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP)         :: var(:)
    real(RP)         :: valmin
    real(RP)         :: valmax
    character(len=*) :: varname

    logical :: invalid_value
    integer :: k, kstr, kend
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug')

    kstr = lbound( var(:), 1 )
    kend = ubound( var(:), 1 )

    invalid_value = .false.
    do k = kstr, kend
       if (      var(k)*0.0_RP /= 0.0_RP &
            .OR. var(k)        <  valmin &
            .OR. var(k)        >  valmax ) then
           invalid_value = .true.
           exit
       endif
    enddo

    if ( invalid_value ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx invalid value:', trim(varname), &
                                      '(', k, ')=', var(k)
       call PRC_MPIstop
    endif

    call PROF_rapend  ('Debug')

    return
  end subroutine MISC_valcheck_1D

  !-----------------------------------------------------------------------------
  !> Nan & extreme value checker (2D)
  subroutine MISC_valcheck_2D( &
       var,    &
       valmin, &
       valmax, &
       varname )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP)         :: var(:,:)
    real(RP)         :: valmin
    real(RP)         :: valmax
    character(len=*) :: varname

    logical :: invalid_value
    integer :: k, kstr, kend
    integer :: i, istr, iend
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug')

    kstr = lbound( var(:,:), 1 )
    kend = ubound( var(:,:), 1 )

    istr = lbound( var(:,:), 2 )
    iend = ubound( var(:,:), 2 )

    invalid_value = .false.
    outer:do i = istr, iend
          do k = kstr, kend
             if (      var(k,i)*0.0_RP /= 0.0_RP &
                  .OR. var(k,i)        <  valmin &
                  .OR. var(k,i)        >  valmax ) then
                 invalid_value = .true.
                 exit outer
             endif
          enddo
          enddo outer

    if ( invalid_value ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx invalid value: ', trim(varname), &
                                      '(', k, ',', i, ')=', var(k,i)
       call PRC_MPIstop
    endif

    call PROF_rapend  ('Debug')

    return
  end subroutine MISC_valcheck_2D

  !-----------------------------------------------------------------------------
  !> Nan & extreme value checker (3D)
  subroutine MISC_valcheck_3D( &
       var,    &
       valmin, &
       valmax, &
       varname )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP)         :: var(:,:,:)
    real(RP)         :: valmin
    real(RP)         :: valmax
    character(len=*) :: varname

    logical :: invalid_value
    integer :: k, kstr, kend
    integer :: i, istr, iend
    integer :: j, jstr, jend
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug')

    kstr = lbound( var(:,:,:), 1 )
    kend = ubound( var(:,:,:), 1 )

    istr = lbound( var(:,:,:), 2 )
    iend = ubound( var(:,:,:), 2 )

    jstr = lbound( var(:,:,:), 3 )
    jend = ubound( var(:,:,:), 3 )

    invalid_value = .false.
    outer:do j = jstr, jend
          do i = istr, iend
          do k = kstr, kend
             if (      var(k,i,j)*0.0_RP /= 0.0_RP &
                  .OR. var(k,i,j)        <  valmin &
                  .OR. var(k,i,j)        >  valmax ) then
                 invalid_value = .true.
                 exit outer
             endif
          enddo
          enddo
          enddo outer

    if ( invalid_value ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx invalid value: ', trim(varname), &
                                      '(', k, ',', i, ',', j, ')=', var(k,i,j)
       call PRC_MPIstop
    endif

    call PROF_rapend  ('Debug')

    return
  end subroutine MISC_valcheck_3D

end module mod_misc
!-------------------------------------------------------------------------------
