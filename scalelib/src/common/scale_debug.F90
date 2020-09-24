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
       KA, KS, KE,   &
       var,          &
       valmin,       &
       valmax,       &
       varname,      &
       current_file, &
       current_line, &
       mask          )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    implicit none
    integer ,         intent(in) :: KA, KS, KE
    real(RP),         intent(in) :: var(KA)
    real(RP),         intent(in) :: valmin
    real(RP),         intent(in) :: valmax
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: current_file
    integer,          intent(in) :: current_line

    logical, intent(in), optional :: mask(KA)

    logical :: invalid_value
    integer :: k
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug', 1)

    invalid_value = .false.
    if ( present(mask) ) then
       do k = KS, KE
          if ( .not. mask(k) ) cycle
          if (      var(k)*0.0_RP /= 0.0_RP &
               .OR. var(k)        <  valmin &
               .OR. var(k)        >  valmax ) then
             invalid_value = .true.
             exit
          endif
       enddo
    else
       do k = KS, KE
          if (      var(k)*0.0_RP /= 0.0_RP &
               .OR. var(k)        <  valmin &
               .OR. var(k)        >  valmax ) then
             invalid_value = .true.
             exit
          endif
       enddo
    end if

    if ( invalid_value ) then
       LOG_ERROR("VALCHECK_1D",*) 'invalid value): ', trim(varname), &
                  '(', k, ')=', var(k)
       LOG_ERROR_CONT(*) 'in file   : ', trim(current_file), ', at line : ', current_line
       call PRC_abort
    endif

    call PROF_rapend  ('Debug', 1)

    return
  end subroutine VALCHECK_1D

  !-----------------------------------------------------------------------------
  !> Nan & extreme value checker (2D)
  subroutine VALCHECK_2D( &
       IA, IS, IE, JA, JS, JE, &
       var,          &
       valmin,       &
       valmax,       &
       varname,      &
       current_file, &
       current_line, &
       mask          )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    implicit none
    integer ,         intent(in) :: IA, IS, IE
    integer ,         intent(in) :: JA, JS, JE
    real(RP),         intent(in) :: var(IA,JA)
    real(RP),         intent(in) :: valmin
    real(RP),         intent(in) :: valmax
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: current_file
    integer,          intent(in) :: current_line

    logical, intent(in), optional :: mask(IA,JA)

    logical :: invalid_value
    integer :: i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug', 1)

    invalid_value = .false.
    if ( present(mask) ) then
       outer1:do j = JS, JE
              do i = IS, IE
                 if ( .not. mask(i,j) ) cycle
                 if (      var(i,j)*0.0_RP /= 0.0_RP &
                      .OR. var(i,j)        <  valmin &
                      .OR. var(i,j)        >  valmax ) then
                    invalid_value = .true.
                    exit outer1
                 endif
              enddo
              enddo outer1
    else
       outer2:do j = JS, JE
              do i = IS, IE
                 if (      var(i,j)*0.0_RP /= 0.0_RP &
                      .OR. var(i,j)        <  valmin &
                      .OR. var(i,j)        >  valmax ) then
                    invalid_value = .true.
                    exit outer2
                 endif
              enddo
              enddo outer2
    end if

    if ( invalid_value ) then
       LOG_ERROR("VALCHECK_2D",*) 'invalid value:', trim(varname), &
                  '(', i, ',', j, ')=', var(i,j)
       LOG_ERROR_CONT(*) 'in file   : ', trim(current_file), ', at line : ', current_line
       call PRC_abort
    endif

    call PROF_rapend  ('Debug', 1)

    return
  end subroutine VALCHECK_2D

  !-----------------------------------------------------------------------------
  !> Nan & extreme value checker (3D)
  subroutine VALCHECK_3D( &
       KA, KS, KE,   &
       IA, IS, IE,   &
       JA, JS, JE,   &
       var,          &
       valmin,       &
       valmax,       &
       varname,      &
       current_file, &
       current_line, &
       mask          )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    implicit none
    integer ,         intent(in) :: KA, KS, KE
    integer ,         intent(in) :: IA, IS, IE
    integer ,         intent(in) :: JA, JS, JE
    real(RP),         intent(in) :: var(KA,IA,JA)
    real(RP),         intent(in) :: valmin
    real(RP),         intent(in) :: valmax
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: current_file
    integer,          intent(in) :: current_line

    logical, intent(in), optional :: mask(IA,JA)

    logical :: invalid_value
    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug', 1)

    invalid_value = .false.
    if ( present(mask) ) then
       outer1:do j = JS, JE
              do i = IS, IE
                 if ( .not. mask(i,j) ) cycle
                 do k = KS, KE
                    if (      var(k,i,j)*0.0_RP /= 0.0_RP &
                         .OR. var(k,i,j)        <  valmin &
                         .OR. var(k,i,j)        >  valmax ) then
                       invalid_value = .true.
                       exit outer1
                    endif
                 enddo
              enddo
              enddo outer1
    else
       outer2:do j = JS, JE
              do i = IS, IE
              do k = KS, KE
                 if (      var(k,i,j)*0.0_RP /= 0.0_RP &
                      .OR. var(k,i,j)        <  valmin &
                      .OR. var(k,i,j)        >  valmax ) then
                    invalid_value = .true.
                    exit outer2
                 endif
              enddo
              enddo
              enddo outer2
    end if

    if ( invalid_value ) then
       LOG_ERROR("VALCHECK_3D",*) 'Invalid value:', trim(varname), &
                  '(', k, ',', i, ',', j, ')=', var(k,i,j)
       LOG_ERROR_CONT(*) 'in file   : ', trim(current_file), ', at line : ', current_line
       call PRC_abort
    endif

    call PROF_rapend  ('Debug', 1)

    return
  end subroutine VALCHECK_3D

end module scale_debug
