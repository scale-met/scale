!-----------------------------------------------------------------
!
!  Module : Spline interpolation for 1-dimensional data
!
!-----------------------------------------------------------------
!
!--- History
!
!    2006/1/10   K.Suzuki  re-constructed for F90 version
!
!-----------------------------------------------------------------
module mod_spline

implicit none
private

public :: getknot, fbspl, fpb, getmatrx, getcoef, fspline
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
  function fpb( ndat, i, kdeg, qknt, xdat, elm, x )

  real :: fpb
  integer :: ndat, i, kdeg
  real(8) :: qknt( ndat+kdeg ), xdat( ndat ), elm( ndat,ndat )
  real(8) :: x

  integer :: l
  real(8) :: sum

  sum = 0.D0
  do l = 1, ndat
    sum = sum + elm( l,i )*fbspl( ndat, l, kdeg, qknt, xdat, x )
  end do

  fpb = sum

  end function fpb
  !---------------------------------------------------------------
  subroutine getmatrx           &
      ( ndat, kdeg, qknt, xdat, & !--- in
        elm                     ) !--- out

  integer, intent(in) :: ndat
  integer, intent(in) :: kdeg
  real(8), intent(in) :: qknt( ndat+kdeg )
  real(8), intent(in) :: xdat( ndat )

  real(8), intent(out) :: elm( ndat,ndat )

  !--- local
  real(8) :: dt
  integer :: iw( 2*ndat ), i, j, inder
  real(8), parameter :: eps = 0.

  do i = 1, ndat
  do j = 1, ndat
    elm( i,j ) = fbspl( ndat, j, kdeg, qknt, xdat, xdat( i ) )
  end do
  end do

  call TINVSS( ndat, elm, dt, eps, ndat, iw, inder )

  return

  end subroutine getmatrx
  !---------------------------------------------------------------
  subroutine getcoef           &
      ( ndat, kdeg, elm, ydat, & !--- in
        coef                   ) !--- out

  integer, intent(in) :: ndat  !  number of data
  integer, intent(in) :: kdeg  !  degree of Spline + 1
  real(8), intent(in) :: elm( ndat,ndat ) ! matrix ( inverse )
  real(8), intent(in) :: ydat( ndat )  ! data of dependent var.

  real(8), intent(out) :: coef( ndat ) ! expansion coefficient

  !--- local
  integer :: i, j
  real(8) :: sum

  do i = 1, ndat
    sum = 0.D0
    do j = 1, ndat
      sum = sum + elm( i,j )*ydat( j )
    end do
    coef( i ) = sum
  end do

  return

  end subroutine getcoef
  !---------------------------------------------------------------
  function fspline ( ndat, kdeg, coef, qknt, xdat, x )

  integer, intent(in) :: ndat
  integer, intent(in) :: kdeg
  real(8), intent(in) :: coef( ndat )
  real(8), intent(in) :: qknt( ndat+kdeg )
  real(8), intent(in) :: xdat( ndat )
  real(8), intent(in) :: x

  real(8) :: fspline

  !--- local
  real(8) :: sum
  integer :: i

  sum = 0.D0
  do i = 1, ndat
    sum = sum + coef( i )*fbspl( ndat, i, kdeg, qknt, xdat, x )
  end do

  fspline = sum

  return

  end function fspline
  !---------------------------------------------------------------
end module mod_spline
!-----------------------------------------------------------------
