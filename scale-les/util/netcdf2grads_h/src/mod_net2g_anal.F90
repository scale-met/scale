!-------------------------------------------------------------------------------------------
!> module NET2G anal
!!
!! @par Description
!!          Analysis module for post-process of scale
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_net2g_anal
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_net2g_vars
  use mod_net2g_error

  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: simple_analysis

  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------------------

  !> simple analysis
  !-----------------------------------------------------------------------------------------
  subroutine simple_analysis( &
      atype,    & ! [in ]
      is,       & ! [in ]
      ie,       & ! [in ]
      js,       & ! [in ]
      je,       & ! [in ]
      isn,      & ! [in ]
      jsn,      & ! [in ]
      nzn,      & ! [in ]
      indata,   & ! [in ]
      outdata   ) ! [out]
    implicit none

    integer,  intent(in)  :: atype
    integer,  intent(in)  :: is, ie, js, je
    integer,  intent(in)  :: isn, jsn, nzn
    real(DP), intent(in)  :: indata(:,:,:,:)
    real(SP), intent(out) :: outdata(:,:)

    integer :: i, j, ni, nj, k
    real(SP) :: work
    !---------------------------------------------------------------------------

    select case( atype )
    case ( a_slice )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          outdata(i,j) = real( indata(ni,nj,1,1) )
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case ( a_max )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          outdata(i,j) = real( maxval(indata(ni,nj,:,1)) )
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case ( a_min )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          outdata(i,j) = real( minval(indata(ni,nj,:,1)) )
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case ( a_sum )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          work = 0.0D0
          do k = 1, nzn
             work = work + real( indata(ni,nj,k,1) )
          enddo
          outdata(i,j) = work
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case ( a_ave )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          work = 0.0D0
          do k = 1, nzn
             work = work + real( indata(ni,nj,k,1) )
          enddo
          outdata(i,j) = work / real(nzn)
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case default
       call err_abort( 0, __LINE__ )
    end select

    return
  end subroutine simple_analysis


end module mod_net2g_anal
