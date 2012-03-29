!-------------------------------------------------------------------------------
!> module Geometrics
!!
!! @par Description
!!          Geometrical convert from plane cartesian coordinate
!!          to earths sphere
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-03-22 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_geometrics
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_FILECHR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GEOMETRICS_setup
  public :: GEOMETRICS_write
  public :: GEOMETRICS_makearea
  public :: GEOMETRICS_makelonlat

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(8), public, save :: GEOMETRICS_lon (1,IA,JA)
  real(8), public, save :: GEOMETRICS_lat (1,IA,JA)

  real(8), public, save :: GEOMETRICS_area(1,IA,JA)
  real(8), public, save :: GEOMETRICS_vol (KA,IA,JA)
  real(8), public, save :: GEOMETRICS_totarea ! total area   (local)
  real(8), public, save :: GEOMETRICS_totvol  ! total volume (local)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(8), private :: GEOMETRICS_startlonlat(2) = (/ 135.2D0, 34.7D0 /)
  real(8), private :: GEOMETRICS_rotation       = 0.D0

  character(len=IO_FILECHR), private :: GEOMETRICS_OUT_BASENAME = ''

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup Geometrics
  !-----------------------------------------------------------------------------
  subroutine GEOMETRICS_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_GEOMETRICS / &
       GEOMETRICS_startlonlat, &
       GEOMETRICS_rotation,    &
       GEOMETRICS_OUT_BASENAME

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[GEOMETRICS]/Categ[GRID]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_GEOMETRICS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_GEOMETRICS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_GEOMETRICS)

    call GEOMETRICS_makelonlat

    call GEOMETRICS_makearea

    if ( GEOMETRICS_OUT_BASENAME /= '' ) then
       call GEOMETRICS_write
    endif

    return
  end subroutine GEOMETRICS_setup

  !-----------------------------------------------------------------------------
  !> Write lon&lat, control area/volume
  !-----------------------------------------------------------------------------
  subroutine GEOMETRICS_write
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use mod_fileio_h, only: &
       FIO_HMID,   &
       FIO_REAL8
    use mod_fileio, only: &
       FIO_output
    implicit none

    real(8) :: sfc(1,IMAX,JMAX)

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output geometrics file ***'

    bname = trim(GEOMETRICS_OUT_BASENAME)
    desc  = 'SCALE3 GEOMETRICS'

    sfc(1,1:IMAX,1:JMAX) = GEOMETRICS_area(1,IS:IE,JS:JE)
    call FIO_output( sfc(:,:,:), bname, desc, '',               &
                     'AREA', 'Control Area', '', 'm2',          &
                     FIO_REAL8, 'ZSFC', 1, 1, 1, NOWSEC, NOWSEC )

    return
  end subroutine GEOMETRICS_write

  !-----------------------------------------------------------------------------
  !> Generate control area/volume
  !-----------------------------------------------------------------------------
  subroutine GEOMETRICS_makearea
    use mod_grid, only: &
       CDZ => GRID_CDZ, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    GEOMETRICS_totarea     = 0.D0
    GEOMETRICS_area(:,:,:) = 0.D0
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       GEOMETRICS_area(1,i,j) = CDX(i) * CDY(j)
       GEOMETRICS_totarea = GEOMETRICS_totarea + GEOMETRICS_area(1,i,j)
    enddo
    enddo
    enddo

    GEOMETRICS_totvol     = 0.D0
    GEOMETRICS_vol(:,:,:) = 0.D0
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       GEOMETRICS_vol(k,i,j) = CDZ(k) * CDX(i) * CDY(j)
       GEOMETRICS_totvol = GEOMETRICS_totvol + GEOMETRICS_vol(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine GEOMETRICS_makearea

  !-----------------------------------------------------------------------------
  !> Generate lon/lat
  !-----------------------------------------------------------------------------
  subroutine GEOMETRICS_makelonlat
    use mod_const, only: &
       PI      => CONST_PI,     &
       ERADIUS => CONST_ERADIUS
    use mod_grid, only: &
       CZ => GRID_CZ, &
       CX => GRID_CX, &
       CY => GRID_CY
    implicit none

    real(8) :: GRID_rotX(1,IA,JA)
    real(8) :: GRID_rotY(1,IA,JA)
    real(8) :: d2r, r2d, theta

    real(8) :: c(2)
    real(8) :: gno(2)
    real(8) :: sph(2)
    real(8) :: rho, gmm

    integer :: i, j
    !---------------------------------------------------------------------------

    d2r   = PI / 180.D0
    r2d   = 180.D0 / PI
    theta = GEOMETRICS_rotation * d2r

    c(1) = GEOMETRICS_startlonlat(1) * d2r
    c(2) = GEOMETRICS_startlonlat(2) * d2r

    do j = JS, JE
    do i = IS, IE
       GRID_rotX(1,i,j) = CY(j) * cos(theta) - CX(i) * sin(theta)
       GRID_rotY(1,i,j) = CY(j) * sin(theta) + CX(i) * cos(theta)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       gno(1) =  GRID_rotX(1,i,j) / ERADIUS ! angle from north-west corner (0,0)[m]
       gno(2) = -GRID_rotY(1,i,j) / ERADIUS ! direction is north->south

       rho = sqrt( gno(1) * gno(1) + gno(2) * gno(2) )
       gmm = atan( rho )

       if ( rho == 0.D0 ) then
          sph(1) = c(1)
          sph(2) = c(2)
       else
          sph(1) = c(1) &
                 + atan( gno(1)*sin(gmm) / ( rho   *cos(c(2))*cos(gmm) &
                                           - gno(2)*sin(c(2))*sin(gmm) ) )
          sph(2) = asin(        sin(c(2))*cos(gmm)       &
                       + gno(2)*cos(c(2))*sin(gmm) / rho )
       endif

       GEOMETRICS_lon(1,i,j) = sph(1) * r2d
       GEOMETRICS_lat(1,i,j) = sph(2) * r2d
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) ' *** Position on the earth (Local)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f9.5,A,f9.5,A,A,f9.5,A,f9.5,A)') &
                                'NW(',GEOMETRICS_lon(1,IS,JS),',',GEOMETRICS_lat(1,IS,JS),')-', &
                                'NE(',GEOMETRICS_lon(1,IS,JE),',',GEOMETRICS_lat(1,IS,JE),')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
                                '            |                       |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f9.5,A,f9.5,A,A,f9.5,A,f9.5,A)') &
                                'SW(',GEOMETRICS_lon(1,IE,JS),',',GEOMETRICS_lat(1,IE,JS),')-', &
                                'SE(',GEOMETRICS_lon(1,IE,JE),',',GEOMETRICS_lat(1,IE,JE),')'

    return
  end subroutine GEOMETRICS_makelonlat

end module mod_geometrics
