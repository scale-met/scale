!-------------------------------------------------------------------------------
!> module vector
!!
!! @par Description
!!          module for 3D vector on the sphere
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_vector
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
  public :: VECTR_xyz2latlon
  public :: VECTR_latlon2xyz
  public :: VECTR_cross
  public :: VECTR_dot
  public :: VECTR_abs
  public :: VECTR_angle
  public :: VECTR_intersec
  public :: VECTR_anticlockwise
  public :: VECTR_triangle
  public :: VECTR_distance

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
  subroutine VECTR_xyz2latlon( &
       x,   &
       y,   &
       z,   &
       lat, &
       lon  )
    use scale_const, only: &
       CONST_EPS
    implicit none

    real(RP),intent(in)  :: x
    real(RP),intent(in)  :: y
    real(RP),intent(in)  :: z
    real(RP),intent(out) :: lat
    real(RP),intent(out) :: lon

    real(RP) :: length, length_h
    !---------------------------------------------------------------------------

    length = sqrt( x*x + y*y + z*z )

    if ( length < CONST_EPS ) then ! 3D vector length is
       lat = 0.0_RP
       lon = 0.0_RP
       return
    endif

    if    ( z / length >= 1.0_RP ) then ! vector is parallele to z axis.
       lat = asin(1.0_RP)
       lon = 0.0_RP
       return
    elseif( z / length <= -1.0_RP ) then ! vector is parallele to z axis.
       lat = asin(-1.0_RP)
       lon = 0.0_RP
       return
    else
       lat = asin( z / length )
    endif

    length_h = sqrt( x*x + y*y )

    if ( length_h < CONST_EPS ) then
       lon = 0.0_RP
       return
    endif

    if    ( x / length_h >= 1.0_RP ) then
       lon = acos(1.0_RP)
    elseif( x / length_h <= -1.0_RP ) then
       lon = acos(-1.0_RP)
    else
       lon = acos( x / length_h )
    endif

    if( y < 0.0_RP ) lon = -lon

    return
  end subroutine VECTR_xyz2latlon

  !-----------------------------------------------------------------------------
  subroutine VECTR_latlon2xyz( &
       lat,   &
       lon,   &
       x,     &
       y,     &
       z,     &
       radius )
    implicit none

    real(RP),intent(in)  :: lat
    real(RP),intent(in)  :: lon
    real(RP),intent(out) :: x
    real(RP),intent(out) :: y
    real(RP),intent(out) :: z
    real(RP),intent(in)  :: radius
    !---------------------------------------------------------------------------

    x = radius * cos(lat) * cos(lon)
    y = radius * cos(lat) * sin(lon)
    z = radius * sin(lat)

    return
  end subroutine VECTR_latlon2xyz

  !-----------------------------------------------------------------------------
  !> exterior product of vector a->b and c->d
  subroutine VECTR_cross( nv, a, b, c, d )
    implicit none

    real(RP), intent(out) :: nv(3)                  ! normal vector
    real(RP), intent(in ) :: a(3), b(3), c(3), d(3) ! x,y,z(cartesian)
    !---------------------------------------------------------------------------

    nv(1) = ( b(2)-a(2) ) * ( d(3)-c(3) ) &
          - ( b(3)-a(3) ) * ( d(2)-c(2) )
    nv(2) = ( b(3)-a(3) ) * ( d(1)-c(1) ) &
          - ( b(1)-a(1) ) * ( d(3)-c(3) )
    nv(3) = ( b(1)-a(1) ) * ( d(2)-c(2) ) &
          - ( b(2)-a(2) ) * ( d(1)-c(1) )

    return
  end subroutine VECTR_cross

  !-----------------------------------------------------------------------------
  !> interior product of vector a->b and c->d
  subroutine VECTR_dot( l, a, b, c, d )
    implicit none

    real(RP), intent(out) :: l
    real(RP), intent(in ) :: a(3), b(3), c(3), d(3) ! x,y,z(cartesian)
    !---------------------------------------------------------------------------
    ! if a=c=zero-vector and b=d, result is abs|a|^2

    l = ( b(1)-a(1) ) * ( d(1)-c(1) ) &
      + ( b(2)-a(2) ) * ( d(2)-c(2) ) &
      + ( b(3)-a(3) ) * ( d(3)-c(3) )

    return
  end subroutine VECTR_dot

  !-----------------------------------------------------------------------------
  !> length of vector o->a
  subroutine VECTR_abs( l, a )
    implicit none

    real(RP), intent(out) :: l
    real(RP), intent(in ) :: a(3) ! x,y,z(cartesian)
    !---------------------------------------------------------------------------

    l = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    l = sqrt(l)

    return
  end subroutine VECTR_abs

  !---------------------------------------------------------------------
  !> calc angle between two vector(b->a,b->c)
  subroutine VECTR_angle( angle, a, b, c )
    implicit none

    real(RP), intent(out) :: angle
    real(RP), intent(in ) :: a(3), b(3), c(3)

    real(RP) :: nv(3), nvlenS, nvlenC
    !---------------------------------------------------------------------

    call VECTR_dot  ( nvlenC, b, a, b, c )
    call VECTR_cross( nv(:),  b, a, b, c )
    call VECTR_abs  ( nvlenS, nv(:) )
    angle = atan2( nvlenS, nvlenC )

    return
  end subroutine VECTR_angle

  !-----------------------------------------------------------------------------
  !> judge intersection of two vector
  subroutine VECTR_intersec( ifcross, p, a, b, c, d )
    implicit none

    logical, intent(out) :: ifcross
    ! .true. : line a->b and c->d intersect
    ! .false.: line a->b and c->d do not intersect and p = (0,0)
    real(RP), intent(out) :: p(3) ! intersection point
    real(RP), intent(in ) :: a(3), b(3), c(3), d(3)

    real(RP), parameter :: o(3) = 0.0_RP

    real(RP)            :: oaob(3), ocod(3), cdab(3)
    real(RP)            :: ip, length
    real(RP)            :: angle_aop, angle_pob, angle_aob
    real(RP)            :: angle_cop, angle_pod, angle_cod

    real(RP), parameter :: eps = 1.E-12_RP
    !---------------------------------------------------------------------

    call VECTR_cross( oaob, o, a, o, b )
    call VECTR_cross( ocod, o, c, o, d )
    call VECTR_cross( cdab, o, ocod, o, oaob )

    call VECTR_abs  ( length, cdab )
    call VECTR_dot  ( ip, o, cdab, o, a )

    p(:) = cdab(:) / sign(length,ip)
!    write(ADM_LOG_FID,*), "p:", p(:)

    call VECTR_angle( angle_aop, a, o, p )
    call VECTR_angle( angle_pob, p, o, b )
    call VECTR_angle( angle_aob, a, o, b )
!    write(ADM_LOG_FID,*), "angle a-p-b:", angle_aop, angle_pob, angle_aob

    call VECTR_angle( angle_cop, c, o, p )
    call VECTR_angle( angle_pod, p, o, d )
    call VECTR_angle( angle_cod, c, o, d )
!    write(ADM_LOG_FID,*), "angle c-p-d:", angle_cop, angle_pod, angle_cod

!    write(ADM_LOG_FID,*), "judge:", angle_aob-(angle_aop+angle_pob), angle_cod-(angle_cop+angle_pod)

    ! --- judge intersection
    if (       abs(angle_aob-(angle_aop+angle_pob)) < eps &
         .AND. abs(angle_cod-(angle_cop+angle_pod)) < eps &
         .AND. abs(angle_aop) > eps                       &
         .AND. abs(angle_pob) > eps                       &
         .AND. abs(angle_cop) > eps                       &
         .AND. abs(angle_pod) > eps                       ) then
       ifcross = .true.
    else
       ifcross = .false.
       p(:) = 0.0_RP
    endif

    return
  end subroutine VECTR_intersec

  !---------------------------------------------------------------------
  !> bubble sort anticlockwise by angle
  subroutine VECTR_anticlockwise( vertex, nvert )
    implicit none

    integer, intent(in)    :: nvert
    real(RP), intent(inout) :: vertex(nvert,3)

    real(RP), parameter :: o(3) = 0.0_RP
    real(RP)            :: v1(3), v2(3), v3(3)
    real(RP)            :: xp(3), ip
    real(RP)            :: angle1, angle2

    real(RP), parameter :: eps = 1.E-12_RP

    integer :: i, j
    !---------------------------------------------------------------------

    do j = 2  , nvert-1
    do i = j+1, nvert
       v1(:) = vertex(1,:)
       v2(:) = vertex(j,:)
       v3(:) = vertex(i,:)

       call VECTR_cross( xp(:), v1(:), v2(:), v1(:), v3(:) )
       call VECTR_dot  ( ip, o(:), v1(:), o(:), xp(:) )

       if ( ip < -eps ) then ! right hand : exchange
!          write(ADM_LOG_FID,*) 'exchange by ip', i, '<->',j
          vertex(i,:) = v2(:)
          vertex(j,:) = v3(:)
       endif

    enddo
    enddo

    v1(:) = vertex(1,:)
    v2(:) = vertex(2,:)
    v3(:) = vertex(3,:)
    ! if 1->2->3 is on the line
    call VECTR_cross( xp(:), v1(:), v2(:), v1(:), v3(:) )
    call VECTR_dot  ( ip, o(:), v1(:), o(:), xp(:) )
    call VECTR_angle( angle1, v1(:), o, v2(:) )
    call VECTR_angle( angle2, v1(:), o, v3(:) )
!    write(ADM_LOG_FID,*) ip, angle1, angle2, abs(angle1)-abs(angle2)

    if (       abs(ip)                 < eps  &      ! on the same line
         .AND. abs(angle2)-abs(angle1) < 0.0_RP ) then ! which is far?
!       write(ADM_LOG_FID,*) 'exchange by angle', 2, '<->', 3
       vertex(2,:) = v3(:)
       vertex(3,:) = v2(:)
    endif

    v2(:) = vertex(nvert  ,:)
    v3(:) = vertex(nvert-1,:)
    ! if 1->nvert->nvert-1 is on the line
    call VECTR_cross( xp(:), v1(:), v2(:), v1(:), v3(:) )
    call VECTR_dot  ( ip, o(:), v1(:), o(:), xp(:) )
    call VECTR_angle( angle1, v1(:), o, v2(:) )
    call VECTR_angle( angle2, v1(:), o, v3(:) )
!    write(ADM_LOG_FID,*) ip, angle1, angle2, abs(angle1)-abs(angle2)

    if (       abs(ip)                 < eps  &      ! on the same line
         .AND. abs(angle2)-abs(angle1) < 0.0_RP ) then ! which is far?
!       write(ADM_LOG_FID,*) 'exchange by angle', nvert, '<->', nvert-1
       vertex(nvert,  :) = v3(:)
       vertex(nvert-1,:) = v2(:)
    endif

    return
  end subroutine VECTR_anticlockwise

  !-----------------------------------------------------------------------------
  function VECTR_triangle( &
       a, b, c,      & !--- IN : three point vectors on a sphere.
       polygon_type, & !--- IN : sphere triangle or plane one?
       radius      ) & !--- IN : radius
       result(area)    !--- OUT : triangle area
    use scale_const, only: &
       CONST_EPS, &
       CONST_PI
    implicit none

    real(RP),         intent(in) :: a(3), b(3), c(3)
    character(len=*), intent(in) :: polygon_type
    real(RP),         intent(in) :: radius
    real(RP)                     :: area

    real(RP), parameter :: o(3) = 0.0_RP

    ! ON_PLANE
    real(RP) :: abc(3)
    real(RP) :: prd, r

    ! ON_SPHERE
    real(RP) :: angle(3)
    real(RP) :: oaob(3), oaoc(3)
    real(RP) :: oboc(3), oboa(3)
    real(RP) :: ocoa(3), ocob(3)
    real(RP) :: abab, acac
    real(RP) :: bcbc, baba
    real(RP) :: caca, cbcb
    !---------------------------------------------------------------------------

    area = 0.0_RP

    if ( polygon_type == 'ON_PLANE' ) then ! Note : On a plane, area = | ourter product of two vectors |.

       call VECTR_cross( abc(:), a(:), b(:), a(:), c(:) )
       call VECTR_abs( prd, abc(:) )
       call VECTR_abs( r  , a(:)   )

       prd = 0.5_RP * prd !! triangle area
       if ( r < CONST_EPS * radius ) then
          print *, "zero length?", a(:)
       else
          r = 1.0_RP / r   !! 1 / length
       endif

       area = prd * r*r * radius*radius

    elseif( polygon_type == 'ON_SPHERE' ) then ! On a unit sphere, area = sum of angles - pi.

       ! angle 1
       call VECTR_cross( oaob(:), o(:), a(:), o(:), b(:) )
       call VECTR_cross( oaoc(:), o(:), a(:), o(:), c(:) )
       call VECTR_abs( abab, oaob(:) )
       call VECTR_abs( acac, oaoc(:) )

       if ( abab < CONST_EPS * radius .OR. acac < CONST_EPS * radius ) then
          !write(*,'(A,3(E20.10))') "zero length abab or acac:", abab/radius, acac/radius
          return
       endif

       call VECTR_angle( angle(1), oaob(:), o(:), oaoc(:) )

       ! angle 2
       call VECTR_cross( oboc(:), o(:), b(:), o(:), c(:) )
       oboa(:) = -oaob(:)
       call VECTR_abs( bcbc, oboc(:) )
       baba = abab

       if ( bcbc < CONST_EPS * radius .OR. baba < CONST_EPS * radius ) then
          !write(*,'(A,3(E20.10))') "zero length bcbc or baba:", bcbc/radius, baba/radius
          return
       endif

       call VECTR_angle( angle(2), oboc(:), o(:), oboa(:) )

       ! angle 3
       ocoa(:) = -oaoc(:)
       ocob(:) = -oboc(:)
       caca = acac
       cbcb = bcbc

       if ( caca < CONST_EPS * radius .OR. cbcb < CONST_EPS * radius ) then
          !write(*,'(A,3(E20.10))') "zero length caca or cbcb:", caca/radius, cbcb/radius
          return
       endif

       call VECTR_angle( angle(3), ocoa(:), o(:), ocob(:) )

       ! calc area
       area = ( angle(1)+angle(2)+angle(3) - CONST_PI ) * radius*radius

    endif

    return
  end function VECTR_triangle

  !-----------------------------------------------------------------------
  !> Get horizontal distance on the sphere
  subroutine VECTR_distance( &
       r,    &
       lon1, &
       lat1, &
       lon2, &
       lat2, &
       dist  )
    implicit none

    real(RP), intent(in)  :: r          ! radius in meter
    real(RP), intent(in)  :: lon1, lat1 ! in radian
    real(RP), intent(in)  :: lon2, lat2 ! in radian
    real(RP), intent(out) :: dist       ! distance of the two points in meter

    real(RP) :: gmm, gno_x, gno_y
    !-----------------------------------------------------------------------

    gmm = sin(lat1) * sin(lat2) &
        + cos(lat1) * cos(lat2) * cos(lon2-lon1)

    gno_x = gmm * ( cos(lat2) * sin(lon2-lon1) )
    gno_y = gmm * ( cos(lat1) * sin(lat2) &
                  - sin(lat1) * cos(lat2) * cos(lon2-lon1) )

    dist = r * atan2( sqrt(gno_x*gno_x+gno_y*gno_y), gmm )

    return
  end subroutine VECTR_distance

end module scale_vector
