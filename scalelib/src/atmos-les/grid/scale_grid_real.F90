!-------------------------------------------------------------------------------
!> module GRID (real space)
!!
!! @par Description
!!          Grid module for orthogonal curvelinear, terrain-following coordinate
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-24 (H.Yashiro)  [new] reconstruction from scale_REAL & scale_topography
!!
!<
!-------------------------------------------------------------------------------
module scale_grid_real
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: REAL_setup
  public :: REAL_update_Z

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: REAL_LON(:,:)       !< longitude [rad,0-2pi]
  real(RP), public, allocatable :: REAL_LAT(:,:)       !< latitude  [rad,-pi,pi]
  real(RP), public, allocatable :: REAL_CZ (:,:,:)     !< geopotential height [m] (cell center)
  real(RP), public, allocatable :: REAL_FZ (:,:,:)     !< geopotential height [m] (cell face  )

  real(RP), public, allocatable :: REAL_LONX(:,:)      !< longitude at staggered point (u) [rad,0-2pi]
  real(RP), public, allocatable :: REAL_LATY(:,:)      !< latitude  at staggered point (v) [rad,-pi,pi]
  real(RP), public, allocatable :: REAL_DLON(:,:)      !< delta longitude
  real(RP), public, allocatable :: REAL_DLAT(:,:)      !< delta latitude

  real(RP), public, allocatable :: REAL_PHI (:,:,:)    !< geopotential [m2/s2] (cell center)

  real(RP), public, allocatable :: REAL_AREA(:,:)      !< horizontal area [m2]
  real(RP), public, allocatable :: REAL_VOL (:,:,:)    !< control volume  [m3]

  real(RP), public :: REAL_TOTAREA                     !< total area   (local) [m2]
  real(RP), public :: REAL_TOTVOL                      !< total volume (local) [m3]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: REAL_calc_latlon
  private :: REAL_calc_Z
  private :: REAL_calc_areavol

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private :: REAL_OUT_BASENAME = ''                     !< basename of the output file
  character(len=H_MID),  private :: REAL_OUT_TITLE    = 'SCALE-LES GEOMETRICS' !< title    of the output file
  character(len=H_MID),  private :: REAL_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine REAL_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_mapproj, only: &
       MPRJ_setup
    use scale_fileio, only: &
       FILEIO_set_coordinates
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[REAL]/Categ[GRID]'

    allocate( REAL_LON (IA,JA) )
    allocate( REAL_LAT (IA,JA) )
    allocate( REAL_LONX(IA,JA) )
    allocate( REAL_LATY(IA,JA) )
    allocate( REAL_DLON(IA,JA) )
    allocate( REAL_DLAT(IA,JA) )

    allocate( REAL_CZ (  KA,IA,JA) )
    allocate( REAL_FZ (0:KA,IA,JA) )
    allocate( REAL_PHI(  KA,IA,JA) )

    allocate( REAL_AREA(IA,JA)    )
    allocate( REAL_VOL (KA,IA,JA) )

    ! setup map projection
    call MPRJ_setup

    ! calc longitude & latitude
    call REAL_calc_latlon

    ! calc real height
    call REAL_calc_Z

    ! calc control area & volume
    call REAL_calc_areavol

    ! set latlon and z to fileio module
    call FILEIO_set_coordinates( REAL_LON, REAL_LONX, &
                                 REAL_LAT, REAL_LATY, &
                                 REAL_CZ,  REAL_FZ    )

    return
  end subroutine REAL_setup

  !-----------------------------------------------------------------------------
  !> Re-setup with updated topography
  subroutine REAL_update_Z
    use scale_process, only: &
       PRC_MPIstop
    use scale_fileio, only: &
       FILEIO_set_coordinates
    implicit none
    !---------------------------------------------------------------------------

    ! calc real height
    call REAL_calc_Z

    ! calc control area & volume
    call REAL_calc_areavol

    ! set latlon and z to fileio module
    call FILEIO_set_coordinates( REAL_LON, REAL_LONX, &
                                 REAL_LAT, REAL_LATY, &
                                 REAL_CZ,  REAL_FZ    )

    return
  end subroutine REAL_update_Z

  !-----------------------------------------------------------------------------
  !> Calc longitude & latitude
  subroutine REAL_calc_latlon
    use scale_const, only: &
       D2R    => CONST_D2R
    use scale_grid, only: &
       CX => GRID_CX, &
       CY => GRID_CY, &
       FX => GRID_FX, &
       FY => GRID_FY
    use scale_mapproj, only: &
       MPRJ_xy2lonlat
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       call MPRJ_xy2lonlat( CX(i), FY(j), REAL_LON(i,j), REAL_LATY(i,j))
       call MPRJ_xy2lonlat( FX(i), CY(j), REAL_LONX(i,j), REAL_LAT(i,j))
    enddo
    enddo

    REAL_DLON(:,:) = 0.0_RP
    REAL_DLAT(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       REAL_DLON(i,j) = REAL_LONX(i,j) - REAL_LONX(i-1,j)
       REAL_DLAT(i,j) = REAL_LATY(i,j) - REAL_LATY(i,j-1)
    enddo
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
  end subroutine REAL_calc_latlon

  !-----------------------------------------------------------------------------
  !> Convert Xi to Z coordinate
  subroutine REAL_calc_Z
    use scale_const, only: &
       CONST_GRAV
    use scale_grid, only: &
       CZ => GRID_CZ, &
       FZ => GRID_FZ
    use scale_topography, only: &
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
  end subroutine REAL_calc_Z

  !-----------------------------------------------------------------------------
  !> Calc control area/volume
  subroutine REAL_calc_areavol
    use scale_const, only: &
       RADIUS => CONST_RADIUS
    use scale_grid, only: &
       DZ, &
       DX, &
       DY
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( .false. ) then
       REAL_TOTAREA   = 0.0_RP
       REAL_AREA(:,:) = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          REAL_AREA(i,j) = RADIUS * RADIUS * REAL_DLON(i,j) &
                         * ( sin( REAL_LAT(i,j)-0.5_RP*REAL_DLAT(i,j) ) &
                           - sin( REAL_LAT(i,j)+0.5_RP*REAL_DLAT(i,j) ) )
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
  end subroutine REAL_calc_areavol

end module scale_grid_real
