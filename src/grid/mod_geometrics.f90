!-------------------------------------------------------------------------------
!> module GEOMETRICS
!!
!! @par Description
!!          Geometrical convert from plane cartesian coordinate
!!          to planet sphere
!!
!! @author Team SCALE
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
     IO_FILECHR, &
     IO_SYSCHR
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include "inc_index.h"

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
  !++ Public parameters & variables
  !
  real(RP), public, save :: GEOMETRICS_lon (1,IA,JA)  !< longitude [rad,0-2pi]
  real(RP), public, save :: GEOMETRICS_lat (1,IA,JA)  !< latitude  [rad,-pi,pi]

  real(RP), public, save :: GEOMETRICS_area(1,IA,JA)  !< cell horizontal area [m2]
  real(RP), public, save :: GEOMETRICS_vol (KA,IA,JA) !< cell volume          [m3]

  real(RP), public, save :: GEOMETRICS_totarea        !< total area   (local) [m2]
  real(RP), public, save :: GEOMETRICS_totvol         !< total volume (local) [m3]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: GEOMETRICS_startlonlat(2) = (/ 135.2_RP, 34.7_RP /) !< lon&lat at north-west corner
  real(RP), private :: GEOMETRICS_rotation       = 0.0_RP                  !< rotation angle

  character(len=IO_FILECHR), private :: GEOMETRICS_OUT_BASENAME = ''                  !< basename of the output file
  character(len=IO_SYSCHR),  private :: GEOMETRICS_OUT_TITLE    = 'SCALE3 GEOMETRICS' !< title    of the output file
  character(len=IO_SYSCHR),  private :: GEOMETRICS_OUT_DTYPE    = 'DEFAULT'           !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
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
       GEOMETRICS_OUT_DTYPE

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

    ! calc metrics
    call GEOMETRICS_makearea
    call GEOMETRICS_makelonlat

    ! write to file
    call GEOMETRICS_write

    return
  end subroutine GEOMETRICS_setup

  !-----------------------------------------------------------------------------
  !> Write lon&lat, control area/volume
  subroutine GEOMETRICS_write
    use mod_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( GEOMETRICS_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output geometrics file ***'

       call FILEIO_write( GEOMETRICS_lon(1,:,:),  GEOMETRICS_OUT_BASENAME, GEOMETRICS_OUT_TITLE, &
                          'lon', 'Longitude', 'degrees_east', 'XY',        GEOMETRICS_OUT_DTYPE  )

       call FILEIO_write( GEOMETRICS_lat(1,:,:),  GEOMETRICS_OUT_BASENAME, GEOMETRICS_OUT_TITLE, &
                          'lat', 'Latitude', 'degrees_north', 'XY',        GEOMETRICS_OUT_DTYPE  )

       call FILEIO_write( GEOMETRICS_area(1,:,:), GEOMETRICS_OUT_BASENAME, GEOMETRICS_OUT_TITLE, &
                          'area', 'Control Area', 'm2', 'XY',              GEOMETRICS_OUT_DTYPE  )

       call FILEIO_write( GEOMETRICS_vol(:,:,:),  GEOMETRICS_OUT_BASENAME, GEOMETRICS_OUT_TITLE, &
                          'vol', 'Control Volume', 'm3', 'ZXY',            GEOMETRICS_OUT_DTYPE  )

    endif

    return
  end subroutine GEOMETRICS_write

  !-----------------------------------------------------------------------------
  !> Calc control area/volume
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
  !> Calc lon/lat
  subroutine GEOMETRICS_makelonlat
    use mod_const, only: &
       RADIUS => CONST_RADIUS, &
       D2R    => CONST_D2R
    use mod_grid, only: &
       CX => GRID_CX, &
       CY => GRID_CY
    implicit none

    real(RP) :: GRID_rotX(1,IA,JA)
    real(RP) :: GRID_rotY(1,IA,JA)
    real(RP) :: theta

    real(RP) :: c(2)
    real(RP) :: gno(2)
    real(RP) :: sph(2)
    real(RP) :: rho, gmm

    integer :: i, j
    !---------------------------------------------------------------------------

    theta = GEOMETRICS_rotation * D2R

    c(1) = GEOMETRICS_startlonlat(1) * D2R
    c(2) = GEOMETRICS_startlonlat(2) * D2R

    do j = JS, JE
    do i = IS, IE
       GRID_rotX(1,i,j) = CY(j) * cos(theta) - CX(i) * sin(theta)
       GRID_rotY(1,i,j) = CY(j) * sin(theta) + CX(i) * cos(theta)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       gno(1) =  GRID_rotX(1,i,j) / RADIUS ! angle from north-west corner (0,0)[m]
       gno(2) = -GRID_rotY(1,i,j) / RADIUS ! direction is north->south

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

       GEOMETRICS_lon(1,i,j) = sph(1)
       GEOMETRICS_lat(1,i,j) = sph(2)
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) ' *** Position on the earth (Local)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f9.5,A,f9.5,A,A,f9.5,A,f9.5,A)') &
                                'NW(',GEOMETRICS_lon(1,IS,JS)/D2R,',',GEOMETRICS_lat(1,IS,JS)/D2R,')-', &
                                'NE(',GEOMETRICS_lon(1,IS,JE)/D2R,',',GEOMETRICS_lat(1,IS,JE)/D2R,')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
                                '            |                       |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f9.5,A,f9.5,A,A,f9.5,A,f9.5,A)') &
                                'SW(',GEOMETRICS_lon(1,IE,JS)/D2R,',',GEOMETRICS_lat(1,IE,JS)/D2R,')-', &
                                'SE(',GEOMETRICS_lon(1,IE,JE)/D2R,',',GEOMETRICS_lat(1,IE,JE)/D2R,')'

    return
  end subroutine GEOMETRICS_makelonlat

end module mod_geometrics
