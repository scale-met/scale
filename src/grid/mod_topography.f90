!-------------------------------------------------------------------------------
!> module TOPOGRAPHY
!!
!! @par Description
!!          Topography module
!!          Terrain-following and vertical metrics
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_topography
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
  public :: TOPO_setup
  public :: TOPO_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: TOPO_Zsfc(IA,JA)    !< absolute ground height [m]

  real(RP), public, save :: TOPO_CZ(  KA,IA,JA) !< Z coordinate [m] (cell center)
  real(RP), public, save :: TOPO_FZ(0:KA,IA,JA) !< Z coordinate [m] (cell face  )

  real(RP), public, save :: TOPO_PHI(KA,IA,JA) !< geopotential (cell center)

  integer,  public, save :: I_XYZ = 1 ! at (x,y,z)
  integer,  public, save :: I_XYW = 2 ! at (x,y,w)
  integer,  public, save :: I_UYW = 3 ! at (u,y,w)
  integer,  public, save :: I_XVW = 4 ! at (x,v,w)
  integer,  public, save :: I_UYZ = 5 ! at (u,y,z)
  integer,  public, save :: I_XVZ = 6 ! at (x,v,z)
  integer,  public, save :: I_UVZ = 7 ! at (u,v,z)

  real(RP), public, save :: TRANSGRID_GSQRT(KA,IA,JA,7) !< transformation metrics from Z to Xi, {G}^1/2
  real(RP), public, save :: TRANSGRID_J13G (KA,IA,JA,4) !< (1,3) element of Jacobian matrix * {G}^1/2
  real(RP), public, save :: TRANSGRID_J23G (KA,IA,JA,4) !< (2,3) element of Jacobian matrix * {G}^1/2
  real(RP), public, save :: TRANSGRID_J33G              !< (3,3) element of Jacobian matrix * {G}^1/2

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TOPO_read
  private :: TOPO_Xi2Z
  private :: TOPO_metrics

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: XDIR = 1 !< [index] X direction
  integer, private, parameter :: YDIR = 2 !< [index] Y direction

  character(len=IO_FILECHR), private :: TOPO_IN_BASENAME  = ''                  !< basename of the input  file
  character(len=IO_FILECHR), private :: TOPO_OUT_BASENAME = ''                  !< basename of the output file
  character(len=IO_SYSCHR),  private :: TOPO_OUT_TITLE    = 'SCALE3 TOPOGRAPHY' !< title    of the output file
  character(len=IO_SYSCHR),  private :: TOPO_OUT_DTYPE    = 'DEFAULT'           !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine TOPO_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_GRAV
    implicit none

    namelist / PARAM_TOPO / &
       TOPO_IN_BASENAME,  &
       TOPO_OUT_BASENAME, &
       TOPO_OUT_DTYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TOPOGRAPHY]/Categ[GRID]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TOPO,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_TOPO. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_TOPO)

    ! read from file
    call TOPO_read

    ! calc metrics
    call TOPO_Xi2Z
    call TOPO_metrics

    TOPO_PHI(:,:,:) = TOPO_CZ(:,:,:) * CONST_GRAV

    ! write to file
    call TOPO_write

    return
  end subroutine TOPO_setup

  !-----------------------------------------------------------------------------
  !> Read topography
  subroutine TOPO_read
    use mod_fileio, only: &
       FILEIO_read
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input topography file ***'

    if ( TOPO_IN_BASENAME /= '' ) then

       call FILEIO_read( TOPO_Zsfc(:,:),                        & ! [OUT]
                         TOPO_IN_BASENAME, 'TOPO', 'XY', step=1 ) ! [IN]
       ! fill IHALO & JHALO
       call COMM_vars8( TOPO_Zsfc(:,:), 1 )
       call COMM_wait ( TOPO_Zsfc(:,:), 1 )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** topography file is not specified.'

       TOPO_Zsfc(:,:) = 0.0_RP
    endif

    return
  end subroutine TOPO_read

  !-----------------------------------------------------------------------------
  !> Write topography
  subroutine TOPO_write
    use mod_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( TOPO_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output topography file ***'

       call FILEIO_write( TOPO_Zsfc(:,:), TOPO_OUT_BASENAME, TOPO_OUT_TITLE, & ! [IN]
                          'TOPO', 'Topography', 'm', 'XY', TOPO_OUT_DTYPE    ) ! [IN]

    endif

    return
  end subroutine TOPO_write

  !-----------------------------------------------------------------------------
  !> Convert Xi to Z coordinate
  subroutine TOPO_Xi2Z
    use mod_grid, only: &
       GRID_CZ, &
       GRID_FZ
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: Htop

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    Htop = GRID_FZ(KE) - GRID_FZ(KS-1)

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       TOPO_CZ(k,i,j) = ( Htop - TOPO_Zsfc(i,j) ) / Htop * GRID_CZ(k) + TOPO_Zsfc(i,j)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 0, KA
       TOPO_FZ(k,i,j) = ( Htop - TOPO_Zsfc(i,j) ) / Htop * GRID_FZ(k) + TOPO_Zsfc(i,j)
    enddo
    enddo
    enddo

    call COMM_vars8( TOPO_CZ(:,:,:), 1 )
    call COMM_vars8( TOPO_FZ(:,:,:), 2 )
    call COMM_wait ( TOPO_CZ(:,:,:), 1 )
    call COMM_wait ( TOPO_FZ(:,:,:), 2 )

    return
  end subroutine TOPO_Xi2Z

  !-----------------------------------------------------------------------------
  !> Calculate G^1/2 & Jacobian
  subroutine TOPO_metrics
    use mod_grid, only: &
       GRID_CX,   &
       GRID_CY,   &
       GRID_CZ,   &
       GRID_FZ,   &
       GRID_RCDZ, &
       GRID_RCDX, &
       GRID_RCDY, &
       GRID_RFDZ, &
       GRID_RFDX, &
       GRID_RFDY
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_fileio, only: &
       FILEIO_write
    implicit none

    real(RP) :: TOPO_CZ_U (  KA,IA,JA) !< Z coordinate [m] at (u,y,z)
    real(RP) :: TOPO_CZ_V (  KA,IA,JA) !< Z coordinate [m] at (x,v,z)
    real(RP) :: TOPO_FZ_U (0:KA,IA,JA) !< Z coordinate [m] at (u,y,w)
    real(RP) :: TOPO_FZ_V (0:KA,IA,JA) !< Z coordinate [m] at (x,v,w)
    real(RP) :: TOPO_FZ_UV(0:KA,IA,JA) !< Z coordinate [m] at (u,v,w)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! calc Z-coordinate height at staggered position
    do j = 1, JA
    do i = 1, IA-1
    do k = 1, KA
       TOPO_CZ_U(k,i,j) = 0.5D0 * ( TOPO_CZ(k,i+1,j) + TOPO_CZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA-1
    do k = 0, KA
       TOPO_FZ_U(k,i,j) = 0.5D0 * ( TOPO_FZ(k,i+1,j) + TOPO_FZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = 1, KA
       TOPO_CZ_V(k,i,j) = 0.5D0 * ( TOPO_CZ(k,i,j+1) + TOPO_CZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = 0, KA
       TOPO_FZ_V(k,i,j) = 0.5D0 * ( TOPO_FZ(k,i,j+1) + TOPO_FZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA-1
    do k = 0, KA
       TOPO_FZ_UV(k,i,j) = 0.25D0 * ( TOPO_FZ(k,i+1,j+1) + TOPO_FZ(k,i+1,j) &
                                    + TOPO_FZ(k,i  ,j+1) + TOPO_FZ(k,i  ,j) )
    enddo
    enddo
    enddo

    ! G^1/2
    do j = JS, JE
    do i = IS, IE
       ! at (x,y,z)
       do k = 1, KA
          TRANSGRID_GSQRT(k,i,j,I_XYZ) = ( TOPO_FZ(k,i,j) - TOPO_FZ(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (x,y,w)
       do k = 1, KA-1
          TRANSGRID_GSQRT(k,i,j,I_XYW) = ( TOPO_CZ(k+1,i,j) - TOPO_CZ(k,i,j) ) * GRID_RFDZ(k)
       enddo

       ! at (u,y,w)
       do k = 1, KA-1
          TRANSGRID_GSQRT(k,i,j,I_UYW) = ( TOPO_CZ_U(k+1,i,j) - TOPO_CZ_U(k,i,j) ) * GRID_RFDZ(k)
       enddo

       ! at (x,v,w)
       do k = 1, KA-1
          TRANSGRID_GSQRT(k,i,j,I_XVW) = ( TOPO_CZ_V(k+1,i,j) - TOPO_CZ_V(k,i,j) ) * GRID_RFDZ(k)
       enddo

       ! at (u,y,z)
       do k = 1, KA
          TRANSGRID_GSQRT(k,i,j,I_UYZ) = ( TOPO_FZ_U(k,i,j) - TOPO_FZ_U(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (x,v,z)
       do k = 1, KA
          TRANSGRID_GSQRT(k,i,j,I_XVZ) = ( TOPO_FZ_V(k,i,j) - TOPO_FZ_V(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (u,v,z)
       do k = 1, KA
          TRANSGRID_GSQRT(k,i,j,I_UVZ) = ( TOPO_FZ_UV(k,i,j) - TOPO_FZ_UV(k-1,i,j) ) * GRID_RCDZ(k)
       enddo
    enddo
    enddo

    call COMM_vars8( TRANSGRID_GSQRT(:,:,:,1), 1 )
    call COMM_vars8( TRANSGRID_GSQRT(:,:,:,2), 2 )
    call COMM_vars8( TRANSGRID_GSQRT(:,:,:,3), 3 )
    call COMM_vars8( TRANSGRID_GSQRT(:,:,:,4), 4 )
    call COMM_vars8( TRANSGRID_GSQRT(:,:,:,5), 5 )
    call COMM_vars8( TRANSGRID_GSQRT(:,:,:,6), 6 )
    call COMM_vars8( TRANSGRID_GSQRT(:,:,:,7), 7 )
    call COMM_wait ( TRANSGRID_GSQRT(:,:,:,1), 1 )
    call COMM_wait ( TRANSGRID_GSQRT(:,:,:,2), 2 )
    call COMM_wait ( TRANSGRID_GSQRT(:,:,:,3), 3 )
    call COMM_wait ( TRANSGRID_GSQRT(:,:,:,4), 4 )
    call COMM_wait ( TRANSGRID_GSQRT(:,:,:,5), 5 )
    call COMM_wait ( TRANSGRID_GSQRT(:,:,:,6), 6 )
    call COMM_wait ( TRANSGRID_GSQRT(:,:,:,7), 7 )

    ! Jacobian * G^1/2
    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       TRANSGRID_J13G(k,i,j,I_XYZ) = -( TOPO_CZ_U (k,i  ,j) - TOPO_CZ_U (k,i-1,j) ) * GRID_RCDX(i)
       TRANSGRID_J13G(k,i,j,I_XYW) = -( TOPO_FZ_U (k,i  ,j) - TOPO_FZ_U (k,i-1,j) ) * GRID_RCDX(i)
       TRANSGRID_J13G(k,i,j,I_UYW) = -( TOPO_FZ   (k,i+1,j) - TOPO_FZ   (k,i  ,j) ) * GRID_RFDX(i)
       TRANSGRID_J13G(k,i,j,I_XVW) = -( TOPO_FZ_UV(k,i  ,j) - TOPO_FZ_UV(k,i-1,j) ) * GRID_RCDX(i)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       TRANSGRID_J23G(k,i,j,I_XYZ) = -( TOPO_CZ_V (k,i,j  ) - TOPO_CZ_V (k,i,j-1) ) * GRID_RCDY(j)
       TRANSGRID_J23G(k,i,j,I_XYW) = -( TOPO_FZ_V (k,i,j  ) - TOPO_FZ_V (k,i,j-1) ) * GRID_RCDY(j)
       TRANSGRID_J23G(k,i,j,I_UYW) = -( TOPO_FZ   (k,i,j+1) - TOPO_FZ   (k,i,j  ) ) * GRID_RFDY(j)
       TRANSGRID_J23G(k,i,j,I_XVW) = -( TOPO_FZ_UV(k,i,j  ) - TOPO_FZ_UV(k,i,j-1) ) * GRID_RCDY(j)
    enddo
    enddo
    enddo

    TRANSGRID_J33G = 1.0_RP ! - 1 / G^1/2 * G^1/2

    ! fill IHALO & JHALO
    call COMM_vars8( TRANSGRID_J13G(:,:,:,I_XYZ),  1 )
    call COMM_vars8( TRANSGRID_J13G(:,:,:,I_XYW),  2 )
    call COMM_vars8( TRANSGRID_J13G(:,:,:,I_UYW),  3 )
    call COMM_vars8( TRANSGRID_J13G(:,:,:,I_XVW),  4 )
    call COMM_vars8( TRANSGRID_J23G(:,:,:,I_XYZ),  5 )
    call COMM_vars8( TRANSGRID_J23G(:,:,:,I_XYW),  6 )
    call COMM_vars8( TRANSGRID_J23G(:,:,:,I_UYW),  7 )
    call COMM_vars8( TRANSGRID_J23G(:,:,:,I_XVW),  8 )

    call COMM_wait ( TRANSGRID_J13G(:,:,:,I_XYZ),  1 )
    call COMM_wait ( TRANSGRID_J13G(:,:,:,I_XYW),  2 )
    call COMM_wait ( TRANSGRID_J13G(:,:,:,I_UYW),  3 )
    call COMM_wait ( TRANSGRID_J13G(:,:,:,I_XVW),  4 )
    call COMM_wait ( TRANSGRID_J23G(:,:,:,I_XYZ),  5 )
    call COMM_wait ( TRANSGRID_J23G(:,:,:,I_XYW),  6 )
    call COMM_wait ( TRANSGRID_J23G(:,:,:,I_UYW),  7 )
    call COMM_wait ( TRANSGRID_J23G(:,:,:,I_XVW),  8 )

    return
  end subroutine TOPO_metrics

end module mod_topography
