!-------------------------------------------------------------------------------
!> module GRID (real space)
!!
!! @par Description
!!          Grid module for orthogonal curvelinear, terrain-following coordinate
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-24 (H.Yashiro)  [new] reconstruction from mod_REAL & mod_topography
!!
!<
!-------------------------------------------------------------------------------
module mod_grid_real
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
  public :: REAL_setup
  public :: REAL_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: REAL_LON(IA,JA)         !< longitude [rad,0-2pi]
  real(RP), public, save :: REAL_LAT(IA,JA)         !< latitude  [rad,-pi,pi]
  real(RP), public, save :: REAL_CZ (  KA,IA,JA)    !< geopotential height [m] (cell center)
  real(RP), public, save :: REAL_FZ (0:KA,IA,JA)    !< geopotential height [m] (cell face  )

!  real(RP), public, save :: REAL_CXYZ(KA,IA,JA,3)   !< absolute position from sphere center [m] (cell center)
!  real(RP), public, save :: REAL_FXYZ(KA,IA,JA,3,3) !< absolute position from sphere center [m] (cell face)

  real(RP), public, save :: REAL_PHI (KA,IA,JA)     !< geopotential [m2/s2] (cell center)

  real(RP), public, save :: REAL_AREA(IA,JA)        !< horizontal area [m2]
  real(RP), public, save :: REAL_VOL (KA,IA,JA)     !< control volume  [m3]

  real(RP), public, save :: REAL_TOTAREA            !< total area   (local) [m2]
  real(RP), public, save :: REAL_TOTVOL             !< total volume (local) [m3]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: REAL_make_latlon
  private :: REAL_make_Z
  private :: REAL_make_area

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), public, save :: REAL_STD_LON = 135.2_RP !< longitude at south-west corner [deg]
  real(RP), public, save :: REAL_STD_LAT =  34.7_RP !< latitude  at south-west corner [deg]
  real(RP), public, save :: REAL_DLON               !< delta longitude
  real(RP), public, save :: REAL_DLAT               !< delta latitude

  character(len=IO_FILECHR), private :: REAL_OUT_BASENAME = ''                  !< basename of the output file
  character(len=IO_SYSCHR),  private :: REAL_OUT_TITLE    = 'SCALE3 GEOMETRICS' !< title    of the output file
  character(len=IO_SYSCHR),  private :: REAL_OUT_DTYPE    = 'DEFAULT'           !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine REAL_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_REAL / &
       REAL_STD_LON,      &
       REAL_STD_LAT,      &
       REAL_OUT_BASENAME, &
       REAL_OUT_DTYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[REAL]/Categ[GRID]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_REAL,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_REAL. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_REAL)

    ! calc longitude & latitude
    call REAL_make_latlon

    ! calc real height
    call REAL_make_Z

    ! calc control area & volume
    call REAL_make_area

    ! write to file
    call REAL_write

    return
  end subroutine REAL_setup

  !-----------------------------------------------------------------------------
  !> Write lon&lat, control area/volume
  subroutine REAL_write
    use mod_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( REAL_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output GEOMETRICAL PARAMETER ***'

       call FILEIO_write( REAL_LON(:,:), REAL_OUT_BASENAME, REAL_OUT_TITLE,        &
                          'lon', 'Longitude', 'degrees_east', 'XY', REAL_OUT_DTYPE )

       call FILEIO_write( REAL_LAT(:,:), REAL_OUT_BASENAME, REAL_OUT_TITLE,        &
                          'lat', 'Latitude', 'degrees_north', 'XY', REAL_OUT_DTYPE )

       call FILEIO_write( REAL_AREA(:,:), REAL_OUT_BASENAME, REAL_OUT_TITLE, &
                          'area', 'Control Area', 'm2', 'XY', REAL_OUT_DTYPE )

       call FILEIO_write( REAL_VOL(:,:,:), REAL_OUT_BASENAME, REAL_OUT_TITLE,  &
                          'vol', 'Control Volume', 'm3', 'ZXY', REAL_OUT_DTYPE )

    endif

    return
  end subroutine REAL_write

  !-----------------------------------------------------------------------------
  !> Calc longitude & latitude
  subroutine REAL_make_latlon
    use mod_const, only: &
       RADIUS => CONST_RADIUS, &
       D2R    => CONST_D2R
    use mod_grid, only: &
       CX => GRID_CX, &
       CY => GRID_CY
    implicit none

    real(RP) :: STD_LON, STD_LAT

    integer :: i, j
    !---------------------------------------------------------------------------

    STD_LON = REAL_STD_LON * D2R
    STD_LAT = REAL_STD_LAT * D2R

    REAL_DLON = DX / RADIUS * cos(REAL_STD_LAT)
    REAL_DLAT = DY / RADIUS

    if( IO_L ) write(IO_FID_LOG,*) ' *** reference point(south-west corner)[lon,lat]=[',REAL_STD_LON,',',REAL_STD_LAT,']'
    if( IO_L ) write(IO_FID_LOG,*) ' *** delta(longitude)[deg] = ', REAL_DLON / D2R
    if( IO_L ) write(IO_FID_LOG,*) ' *** delta(latitude )[deg] = ', REAL_DLAT / D2R

    do i = 1, IA
       REAL_LON(i,:) = STD_LON + CX(i) / RADIUS * cos(REAL_STD_LAT)
    enddo

    do j = 1, JA
       REAL_LAT(:,j) = STD_LAT + CY(j) / RADIUS
    enddo

    if( IO_L ) write(IO_FID_LOG,*) ' *** Position on the earth (Local)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f9.5,A,f9.5,A,A,f9.5,A,f9.5,A)') &
                                'NW(',REAL_LON(IS,JS)/D2R,',',REAL_LAT(IS,JS)/D2R,')-', &
                                'NE(',REAL_LON(IS,JE)/D2R,',',REAL_LAT(IS,JE)/D2R,')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
                                '            |                       |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f9.5,A,f9.5,A,A,f9.5,A,f9.5,A)') &
                                'SW(',REAL_LON(IE,JS)/D2R,',',REAL_LAT(IE,JS)/D2R,')-', &
                                'SE(',REAL_LON(IE,JE)/D2R,',',REAL_LAT(IE,JE)/D2R,')'

    return
  end subroutine REAL_make_latlon

  !-----------------------------------------------------------------------------
  !> Convert Xi to Z coordinate
  subroutine REAL_make_Z
    use mod_const, only : &
       CONST_GRAV
    use mod_grid, only: &
       CZ => GRID_CZ, &
       FZ => GRID_FZ
    use mod_topography, only: &
       Zsfc => TOPO_Zsfc
    implicit none

    real(RP) :: Htop

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    Htop = FZ(KE) - FZ(KS-1)

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       REAL_CZ(k,i,j) = ( Htop - Zsfc(i,j) ) / Htop * CZ(k) + Zsfc(i,j)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 0, KA
       REAL_FZ(k,i,j) = ( Htop - Zsfc(i,j) ) / Htop * FZ(k) + Zsfc(i,j)
    enddo
    enddo
    enddo

    REAL_PHI(:,:,:) = REAL_CZ(:,:,:) * CONST_GRAV

    return
  end subroutine REAL_make_Z

  !-----------------------------------------------------------------------------
  !> Calc control area/volume
  subroutine REAL_make_area
    use mod_const, only: &
       RADIUS => CONST_RADIUS
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( .false. ) then
       REAL_TOTAREA   = 0.0_RP
       REAL_AREA(:,:) = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          REAL_AREA(i,j) = RADIUS * RADIUS * REAL_DLON &
                         * ( sin( REAL_LAT(i,j)-0.5_RP*REAL_DLAT ) &
                           - sin( REAL_LAT(i,j)+0.5_RP*REAL_DLAT ) )
          REAL_TOTAREA = REAL_TOTAREA + REAL_AREA(i,j)
       enddo
       enddo
       enddo
    else
       REAL_TOTAREA   = 0.0_RP
       REAL_AREA(:,:) = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          REAL_AREA(i,j) = DX * DY
          REAL_TOTAREA = REAL_TOTAREA + REAL_AREA(i,j)
       enddo
       enddo
       enddo
    endif

    REAL_TOTVOL     = 0.0_RP
    REAL_VOL(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       REAL_VOL(k,i,j) = ( REAL_FZ(k,i,j) - REAL_FZ(k-1,i,j) ) * REAL_AREA(i,j)
       REAL_TOTVOL = REAL_TOTVOL + REAL_VOL(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine REAL_make_area

end module mod_grid_real
