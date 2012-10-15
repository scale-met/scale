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
  include "inc_precision.h"
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: GEOMETRICS_lon (1,IA,JA)
  real(RP), public, save :: GEOMETRICS_lat (1,IA,JA)

  real(RP), public, save :: GEOMETRICS_area(1,IA,JA)
  real(RP), public, save :: GEOMETRICS_vol (KA,IA,JA)
  real(RP), public, save :: GEOMETRICS_totarea ! total area   (local)
  real(RP), public, save :: GEOMETRICS_totvol  ! total volume (local)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: GEOMETRICS_startlonlat(2) = (/ 135.2E0_RP, 34.7E0_RP /)
  real(RP), private :: GEOMETRICS_rotation       = 0.E0_RP

  character(len=IO_FILECHR), private :: GEOMETRICS_OUT_BASENAME = ''
  character(len=IO_FILECHR), private :: GEOMETRICS_OUT_TITLE = 'SCALE3 GEOMETRICS'
  character(len=IO_FILECHR), private :: GEOMETRICS_OUT_SOURCE = 'SCALE-LES ver. 3' 
  character(len=IO_FILECHR), private :: GEOMETRICS_OUT_INSTITUTE = 'AICS/RIKEN'

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
       GEOMETRICS_startlonlat,  &
       GEOMETRICS_rotation,     &
       GEOMETRICS_OUT_BASENAME, &
       GEOMETRICS_OUT_TITLE,    &
       GEOMETRICS_OUT_SOURCE,   &
       GEOMETRICS_OUT_INSTITUTE

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
    use mod_process, only: &
       PRC_master, &
       PRC_myrank
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use gtool_file_h, only: &
       File_HMID,  &
       File_REAL4, &
       File_REAL8
    use gtool_file, only: &
       FileCreate, &
       FileAddVariable, &
       FilePutAxis, &
       FileWrite, &
       FileClose
    use mod_grid, only: &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY
    implicit none

    real(RP) :: sfc(1,IMAX,JMAX)
    real(RP) :: var(KMAX,IMAX,JMAX)

    character(len=IO_FILECHR) :: bname
    integer :: fid, vid(4)
    integer, parameter :: I_LON  = 1
    integer, parameter :: I_LAT  = 2
    integer, parameter :: I_AREA = 3
    integer, parameter :: I_VOL  = 4
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output geometrics file ***'

    bname = trim(GEOMETRICS_OUT_BASENAME)
    call FileCreate( fid,                                       & ! (out)
         bname,                                                 & ! (in)
         GEOMETRICS_OUT_TITLE,                                  & ! (in)
         GEOMETRICS_OUT_SOURCE,                                 & ! (in)
         GEOMETRICS_OUT_INSTITUTE,                              & ! (in)
         (/'z','x','y'/), (/KMAX,IMAX,JMAX/), (/'Z','X','Y'/),  & ! (in)
         (/'m','m','m'/), (/File_REAL4,File_REAL4,File_REAL4/), & ! (in)
         PRC_master, PRC_myrank                                 ) ! (in)

    call FileAddVariable( vid(I_LON),             & ! (out)
         fid, 'lon', 'Longitude', 'degrees_east', & ! (in)
         (/'x','y'/), File_REAL8                  ) ! (in)
    call FileAddVariable( vid(I_LAT),             & ! (out)
         fid, 'lat', 'Latitude', 'degrees_north', & ! (in)
         (/'x','y'/), File_REAL8                  ) ! (in)
    call FileAddVariable( vid(I_AREA),            & ! (out)
         fid, 'area', 'Control Area', 'm2',       & ! (in)
         (/'x','y'/), File_REAL8                  ) ! (in)
    call FileAddVariable( vid(I_VOL),             & ! (out)
         fid, 'vol', 'Control Volume', 'm3',      & ! (in)
         (/'z','x','y'/), File_REAL8              ) ! (in)

    call FilePutAxis(fid, 'z', GRID_CZ(KS:KE))
    call FilePutAxis(fid, 'x', GRID_CX(KS:KE))
    call FilePutAxis(fid, 'y', GRID_CY(KS:KE))

    sfc(1,1:IMAX,1:JMAX) = GEOMETRICS_lon(1,IS:IE,JS:JE)
    call FileWrite( vid(I_LON), sfc(1,:,:), NOWSEC, NOWSEC )

    sfc(1,1:IMAX,1:JMAX) = GEOMETRICS_lat(1,IS:IE,JS:JE)
    call FileWrite( vid(I_LAT), sfc(1,:,:), NOWSEC, NOWSEC )

    sfc(1,1:IMAX,1:JMAX) = GEOMETRICS_area(1,IS:IE,JS:JE)
    call FileWrite( vid(I_AREA), sfc(1,:,:), NOWSEC, NOWSEC )

    var(1:KMAX,1:IMAX,1:JMAX) = GEOMETRICS_vol(KS:KE,IS:IE,JS:JE)
    call FileWrite( vid(I_VOL), sfc(:,:,:), NOWSEC, NOWSEC )

    call FileClose( fid )

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

    GEOMETRICS_totarea     = 0.0_RP
    GEOMETRICS_area(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       GEOMETRICS_area(1,i,j) = CDX(i) * CDY(j)
       GEOMETRICS_totarea = GEOMETRICS_totarea + GEOMETRICS_area(1,i,j)
    enddo
    enddo
    enddo

    GEOMETRICS_totvol     = 0.0_RP
    GEOMETRICS_vol(:,:,:) = 0.0_RP
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
       CX => GRID_CX, &
       CY => GRID_CY
    implicit none

    real(RP) :: GRID_rotX(1,IA,JA)
    real(RP) :: GRID_rotY(1,IA,JA)
    real(RP) :: d2r, r2d, theta

    real(RP) :: c(2)
    real(RP) :: gno(2)
    real(RP) :: sph(2)
    real(RP) :: rho, gmm

    integer :: i, j
    !---------------------------------------------------------------------------

    d2r   = PI / 180.0_RP
    r2d   = 180.0_RP / PI
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

       if ( rho == 0.0_RP ) then
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
