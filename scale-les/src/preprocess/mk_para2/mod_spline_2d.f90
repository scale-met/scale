!-----------------------------------------------------------------
!
!  Module : Spline interpolation for 2-dimensional data
!
!-----------------------------------------------------------------
!
!--- History
!
!    2006/1/10   K.Suzuki  re-constructed for F90 version
!    2006/4/6    K.Suzuki  extended into 2D version
!
!-----------------------------------------------------------------
module mod_spline_2d

implicit none
private

public :: getknot, fbspl, getcoef2, fspline2
!-----------------------------------------------------------------
contains
  !---------------------------------------------------------------
  subroutine getknot      &
      ( ndat, kdeg, xdat, & !--- in
        qknt              ) !--- out

  integer, intent(in) :: ndat  !  number of data
  integer, intent(in) :: kdeg  !  degree of Spline + 1
  real(8), intent(in) :: xdat( ndat ) ! data of independent var.

  real(8), intent(out) :: qknt( ndat+kdeg ) ! knots for B-Spline

  !--- local
  integer :: i

  do i = 1, kdeg
    qknt( i ) = xdat( 1 )
  end do

  do i = 1, ndat-kdeg
    qknt( i+kdeg ) = ( xdat( i )+xdat( i+kdeg ) )*0.5D0
  end do

  do i = 1, kdeg
    qknt( ndat+i ) = xdat( ndat )
  end do

  return

  end subroutine getknot
  !---------------------------------------------------------------
  recursive function fbspl ( ndat, inum, kdeg, qknt, xdat, x ) &
  result (bspl)

  real(8) :: bspl

  integer, intent(in) :: ndat  !  number of data
  integer, intent(in) :: inum  !  index of B-Spline
  integer, intent(in) :: kdeg  !  degree of B-Spline + 1
  real(8), intent(in) :: qknt( ndat+kdeg )  !  knot of B-Spline
  real(8), intent(in) :: xdat( ndat ) ! data of independent variable
  real(8), intent(in) :: x     !  interpolation point

  !--- local
  real(8) :: bsp1, bsp2

  if ( ( inum == 1 .and. x == xdat( 1 ) ) .or. &
       ( inum == ndat .and. x == xdat( ndat ) ) ) then
    bspl = 1.
    return
  end if

  if ( kdeg == 1 ) then
    if ( x >= qknt( inum ) .and. x < qknt( inum+1 ) ) then
      bspl = 1.D0
    else
      bspl = 0.D0
    end if
  else
    if ( qknt( inum+kdeg-1 ) /= qknt( inum ) ) then
      bsp1 = ( x-qknt( inum ) ) &
            /( qknt( inum+kdeg-1 )-qknt( inum ) ) &
           * fbspl( ndat, inum, kdeg-1, qknt, xdat, x )
    else
      bsp1 = 0.D0
    end if
    if ( qknt( inum+kdeg ) /= qknt( inum+1 ) ) then
      bsp2 = ( qknt( inum+kdeg )-x ) &
            /( qknt( inum+kdeg )-qknt( inum+1 ) ) &
           * fbspl( ndat, inum+1, kdeg-1, qknt, xdat, x )
    else
      bsp2 =  0.D0
    end if
    bspl = bsp1 + bsp2
  end if

  end function fbspl
  !---------------------------------------------------------------
  subroutine getcoef2                 &
      ( mdat, ndat, kdeg, ldeg,       & !--- in
        xdat, ydat, qknt, rknt, zdat, & !--- in
        coef                          ) !--- out

  integer, intent(in) :: mdat  !  number of data (x-direction)
  integer, intent(in) :: ndat  !  number of data (y-direction)
  integer, intent(in) :: kdeg  !  degree of Spline + 1 (x)
  integer, intent(in) :: ldeg  !  degree of Spline + 1 (y)
  real(8), intent(in) :: xdat( mdat )  ! data of independent var. (x)
  real(8), intent(in) :: ydat( ndat )  ! data of independent var. (y)
  real(8), intent(in) :: qknt( mdat+kdeg ) ! knots of B-Spline (x)
  real(8), intent(in) :: rknt( ndat+ldeg ) ! knots of B-Spline (y)
  real(8), intent(in) :: zdat( mdat,ndat ) ! data of dependent var.

  real(8), intent(out) :: coef( mdat,ndat ) ! expansion coefficient

  !--- local
  real(8) :: elmx( mdat,mdat ), elmy( ndat,ndat )
  integer :: iw1( 2*mdat ), iw2( 2*ndat )
  real(8) :: beta( mdat,ndat ), sum, dt
  real(8), parameter :: eps = 0.D0
  integer :: i, j, k, l, inder


  do i = 1, mdat
  do j = 1, mdat
    elmx( i,j ) = fbspl( mdat, j, kdeg, qknt, xdat, xdat( i ) )
  end do
  end do
  call TINVSS( mdat, elmx, dt, eps, mdat, iw1, inder )

  do l = 1, ndat
  do i = 1, mdat
    sum = 0.D0
    do j = 1, mdat
      sum = sum + elmx( i,j )*zdat( j,l )
    end do
    beta( i,l ) = sum
  end do
  end do

  do i = 1, ndat
  do j = 1, ndat
    elmy( i,j ) = fbspl( ndat, j, ldeg, rknt, ydat, ydat( i ) )
  end do
  end do
  call TINVSS( ndat, elmy, dt, eps, ndat, iw2, inder )

  do k = 1, mdat
  do i = 1, ndat
    sum = 0.D0
    do j = 1, ndat
      sum = sum + elmy( i,j )*beta( k,j )
    end do
    coef( k,i ) = sum
  end do
  end do

  return

  end subroutine getcoef2
  !---------------------------------------------------------------
  function fspline2                          &
             ( mdat, ndat, kdeg, ldeg,       & !--- in
               coef, qknt, rknt, xdat, ydat, & !--- in
               x, y                          ) !--- in

  integer, intent(in) :: mdat
  integer, intent(in) :: ndat
  integer, intent(in) :: kdeg
  integer, intent(in) :: ldeg
  real(8), intent(in) :: coef( mdat,ndat )
  real(8), intent(in) :: qknt( mdat+kdeg )
  real(8), intent(in) :: rknt( ndat+ldeg )
  real(8), intent(in) :: xdat( mdat ), ydat( ndat )
  real(8), intent(in) :: x, y

  real(8) :: fspline2

  !--- local
  real(8) :: sum, add
  integer :: i, j

  sum = 0.D0
  do i = 1, mdat
  do j = 1, ndat
    add = coef( i,j )*fbspl( mdat, i, kdeg, qknt, xdat, x ) &
                     *fbspl( ndat, j, ldeg, rknt, ydat, y )
    sum = sum + add
  end do
  end do

  fspline2 = sum

  return

  end function fspline2
  !---------------------------------------------------------------
end module mod_spline_2d
!-----------------------------------------------------------------
