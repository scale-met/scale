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
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'

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

    integer :: k, kstr, kend
    !---------------------------------------------------------------------------

    kstr = lbound( var(:), 1 )
    kend = ubound( var(:), 1 ) 

    do k = kstr, kend
       if (      var(k)*0.0_RP /= 0.0_RP &
            .OR. var(k)        <  valmin &
            .OR. var(k)        >  valmax ) then

          if( IO_L ) write(IO_FID_LOG,*) 'xxx invalid value:', trim(varname), &
                                         '(', k, ')=', var(k)
          call PRC_MPIstop

       endif
    enddo

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

    integer :: k, kstr, kend
    integer :: i, istr, iend
    !---------------------------------------------------------------------------

    kstr = lbound( var(:,:), 1 )
    kend = ubound( var(:,:), 1 ) 

    istr = lbound( var(:,:), 2 )
    iend = ubound( var(:,:), 2 ) 

    do k = kstr, kend
    do i = istr, iend
       if (      var(k,i)*0.0_RP /= 0.0_RP &
            .OR. var(k,i)        <  valmin &
            .OR. var(k,i)        >  valmax ) then

          if( IO_L ) write(IO_FID_LOG,*) 'xxx invalid value: ', trim(varname), &
                                         '(', k, ',', i, ')=', var(k,i)
          call PRC_MPIstop

       endif
    enddo
    enddo

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

    integer :: k, kstr, kend
    integer :: i, istr, iend
    integer :: j, jstr, jend
    !---------------------------------------------------------------------------

    kstr = lbound( var(:,:,:), 1 )
    kend = ubound( var(:,:,:), 1 ) 

    istr = lbound( var(:,:,:), 2 )
    iend = ubound( var(:,:,:), 2 ) 

    jstr = lbound( var(:,:,:), 3 )
    jend = ubound( var(:,:,:), 3 ) 

    do k = kstr, kend
    do i = istr, iend
    do j = jstr, jend
       if (      var(k,i,j)*0.0_RP /= 0.0_RP &
            .OR. var(k,i,j)        <  valmin &
            .OR. var(k,i,j)        >  valmax ) then

          if( IO_L ) write(IO_FID_LOG,*) 'xxx invalid value: ', trim(varname), &
                                         '(', k, ',', i, ',', j, ')=', var(k,i,j)
          call PRC_MPIstop

       endif
    enddo
    enddo
    enddo

    return
  end subroutine MISC_valcheck_3D

end module mod_misc
!-------------------------------------------------------------------------------
