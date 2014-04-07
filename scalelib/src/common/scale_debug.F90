!-------------------------------------------------------------------------------
!> module DEBUG
!!
!! @par Description
!!          Debug & Value check tools
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-01-15 (H.Yashiro) Merge scale_debug & scale_misc(valcheck)
!!
!<
module scale_debug
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
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
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer,  intent(in) :: current_line
    real(RP), intent(in) :: v
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug')

    if ( .NOT. ( abs(v) < abs(CONST_UNDEF) ) ) then
!       if( IO_L ) write(IO_FID_LOG,,*) 'xxx uninitialized value at line:', current_line

       write(*,*) 'xxx uninitialized value at line:', current_line

       call abort
    end if

    call PROF_rapend  ('Debug')

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
    use scale_process, only: &
       PRC_MPIstop
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
!       if( IO_L ) write(IO_FID_LOG,,*) 'xxx invalid value:', trim(varname), '(', k, ')=', var(k)
!       if( IO_L ) write(IO_FID_LOG,,*) 'xxx in file:', trim(current_file), ', at line:', current_line

       write(*,*) 'xxx invalid value:', trim(varname), '(', k, ')=', var(k)
       write(*,*) 'xxx in file:', trim(current_file), ', at line:', current_line

       call PRC_MPIstop
    endif

    call PROF_rapend  ('Debug')

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
    use scale_process, only: &
       PRC_MPIstop
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
!       if( IO_L ) write(IO_FID_LOG,,*) 'xxx invalid value:', trim(varname), '(', k, ',', i, ')=', var(k,i)
!       if( IO_L ) write(IO_FID_LOG,,*) 'xxx in file:', trim(current_file), ', at line:', current_line

       write(*,*) 'xxx invalid value:', trim(varname), '(', k, ',', i, ')=', var(k,i)
       write(*,*) 'xxx in file:', trim(current_file), ', at line:', current_line

       call PRC_MPIstop
    endif

    call PROF_rapend  ('Debug')

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
    use scale_process, only: &
       PRC_MPIstop
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
!       if( IO_L ) write(IO_FID_LOG,,*) 'xxx invalid value:', trim(varname), '(', k, ',', i, ',', j, ')=', var(k,i,j)
!       if( IO_L ) write(IO_FID_LOG,,*) 'xxx in file:', trim(current_file), ', at line:', current_line

       write(*,*) 'xxx invalid value:', trim(varname), '(', k, ',', i, ',', j, ')=', var(k,i,j)
       write(*,*) 'xxx in file:', trim(current_file), ', at line:', current_line

       call PRC_MPIstop
    endif

    call PROF_rapend  ('Debug')

    return
  end subroutine VALCHECK_3D

end module scale_debug
!-------------------------------------------------------------------------------
