!-------------------------------------------------------------------------------
!> module DEBUG
!!
!! @par Description
!!          Debug & Value check tools
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_debug
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CHECK
  public :: VALCHECK

  interface VALCHECK
     module procedure VALCHECK_1D
     module procedure VALCHECK_2D
     module procedure VALCHECK_3D
  end interface VALCHECK

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  integer, public :: DEBUG_DOMAIN_NUM = 0

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
  !> Undefined value checker
  subroutine CHECK( &
       current_line, &
       v             )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer,  intent(in) :: current_line
    real(RP), intent(in) :: v
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug', 1)

    if ( .NOT. ( abs(v) < abs(CONST_UNDEF) ) ) then
       LOG_ERROR("CHECK",*) 'uninitialized value at line:', current_line
       call abort
    end if

    call PROF_rapend  ('Debug', 1)

    return
  end subroutine CHECK

  !-----------------------------------------------------------------------------
  !> Nan & extreme value checker (1D)
  subroutine VALCHECK_1D( &
       var,          &
       valmin,       &
       valmax,       &
       varname,      &
       current_file, &
       current_line  )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    implicit none

    real(RP),         intent(in) :: var(:)
    real(RP),         intent(in) :: valmin
    real(RP),         intent(in) :: valmax
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: current_file
    integer,          intent(in) :: current_line

    logical :: invalid_value
    integer :: k, kstr, kend
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug', 1)

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
       LOG_ERROR("VALCHECK_1D",*) '[VALCHECK_1D] invalid value:', trim(varname), &
                  '(', PRC_myrank, ',', k, ')=', var(k)
       LOG_ERROR_CONT(*) 'in file   : ', trim(current_file), ', at line : ', current_line
       LOG_ERROR_CONT(*) 'in domain : ', DEBUG_DOMAIN_NUM
       call PRC_abort
    endif

    call PROF_rapend  ('Debug', 1)

    return
  end subroutine VALCHECK_1D

  !-----------------------------------------------------------------------------
  !> Nan & extreme value checker (2D)
  subroutine VALCHECK_2D( &
       var,          &
       valmin,       &
       valmax,       &
       varname,      &
       current_file, &
       current_line  )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    implicit none

    real(RP),         intent(in) :: var(:,:)
    real(RP),         intent(in) :: valmin
    real(RP),         intent(in) :: valmax
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: current_file
    integer,          intent(in) :: current_line

    logical :: invalid_value
    integer :: k, kstr, kend
    integer :: i, istr, iend
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug', 1)

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
       LOG_ERROR("VALCHECK_2D",*) '[VALCHECK_2D] invalid value:', trim(varname), &
                  '(', PRC_myrank, ',', k, ',', i, ')=', var(k,i)
       LOG_ERROR_CONT(*) 'in file   : ', trim(current_file), ', at line : ', current_line
       LOG_ERROR_CONT(*) 'in domain : ', DEBUG_DOMAIN_NUM
       call PRC_abort
    endif

    call PROF_rapend  ('Debug', 1)

    return
  end subroutine VALCHECK_2D

  !-----------------------------------------------------------------------------
  !> Nan & extreme value checker (3D)
  subroutine VALCHECK_3D( &
       var,          &
       valmin,       &
       valmax,       &
       varname,      &
       current_file, &
       current_line  )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    implicit none

    real(RP),         intent(in) :: var(:,:,:)
    real(RP),         intent(in) :: valmin
    real(RP),         intent(in) :: valmax
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: current_file
    integer,          intent(in) :: current_line

    logical :: invalid_value
    integer :: k, kstr, kend
    integer :: i, istr, iend
    integer :: j, jstr, jend
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug', 1)

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
       LOG_ERROR("VALCHECK_3D",*) '[VALCHECK_3D] Invalid value:', trim(varname), &
                  '(', PRC_myrank, ',', k, ',', i, ',', j, ')=', var(k,i,j)
       LOG_ERROR_CONT(*) 'in file   : ', trim(current_file), ', at line : ', current_line
       LOG_ERROR_CONT(*) 'in domain : ', DEBUG_DOMAIN_NUM
       call PRC_abort
    endif

    call PROF_rapend  ('Debug', 1)

    return
  end subroutine VALCHECK_3D

end module scale_debug
