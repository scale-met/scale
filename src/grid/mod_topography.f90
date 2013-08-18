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
  real(RP), public, save :: TOPO_Zsfc  (1,IA,JA)    !< absolute ground height [m]

  real(RP), public, save :: TOPO_CZ    (KA,IA,JA)   !< Z coordinate [m] (cell center)
  real(RP), public, save :: TOPO_FZ    (KA,IA,JA)   !< Z coordinate [m] (cell face  )

  real(RP), public, save :: TOPO_GSQRT (KA,IA,JA)   !< vertical metrics {G}^1/2 at (x,y,layer)
  real(RP), public, save :: TOPO_GSQRTX(KA,IA,JA)   !< vertical metrics {G}^1/2 at (u,y,layer)
  real(RP), public, save :: TOPO_GSQRTY(KA,IA,JA)   !< vertical metrics {G}^1/2 at (x,v,layer)
  real(RP), public, save :: TOPO_GSQRTZ(KA,IA,JA)   !< vertical metrics {G}^1/2 at (x,y,interface)

  real(RP), public, save :: TOPO_J13   (KA,IA,JA)   !< (1,3) element of inverse Jacobian matrix at (x,y,interface)
  real(RP), public, save :: TOPO_J13X  (KA,IA,JA)   !< (1,3) element of inverse Jacobian matrix at (u,y,interface)
  real(RP), public, save :: TOPO_J13Y  (KA,IA,JA)   !< (1,3) element of inverse Jacobian matrix at (x,v,interface)
  real(RP), public, save :: TOPO_J23   (KA,IA,JA)   !< (2,3) element of inverse Jacobian matrix at (x,y,interface)
  real(RP), public, save :: TOPO_J23X  (KA,IA,JA)   !< (2,3) element of inverse Jacobian matrix at (u,y,interface)
  real(RP), public, save :: TOPO_J23Y  (KA,IA,JA)   !< (2,3) element of inverse Jacobian matrix at (x,v,interface)

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

       call FILEIO_read( TOPO_Zsfc(1,:,:),                      & ! [OUT]
                         TOPO_IN_BASENAME, 'TOPO', 'XY', step=1 ) ! [IN]
       ! fill IHALO & JHALO
       call COMM_vars8( TOPO_Zsfc(1,:,:), 1 )
       call COMM_wait ( TOPO_Zsfc(1,:,:), 1 )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** topography file is not specified.'

       TOPO_Zsfc(:,:,:) = 0.0_RP
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

       call FILEIO_write( TOPO_Zsfc(1,:,:),  TOPO_OUT_BASENAME, TOPO_OUT_TITLE, & ! [IN]
                          'TOPO', 'Topography', 'm', 'XY',      TOPO_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine TOPO_write

  !-----------------------------------------------------------------------------
  !> Convert Xi to Z coordinate
  subroutine TOPO_Xi2Z
    use mod_grid, only: &
       GRID_CZ, &
       GRID_FZ
    implicit none

    real(RP) :: Htop

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    Htop = GRID_FZ(KE) - GRID_FZ(KS-1)

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       TOPO_CZ(k,i,j) = ( Htop - TOPO_Zsfc(1,i,j) ) / Htop * GRID_CZ(k) + TOPO_Zsfc(1,i,j)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       TOPO_FZ(k,i,j) = ( Htop - TOPO_Zsfc(1,i,j) ) / Htop * GRID_FZ(k) + TOPO_Zsfc(1,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine TOPO_Xi2Z

  !-----------------------------------------------------------------------------
  !> Calculate G^1/2 & Gvector
  subroutine TOPO_metrics
    use mod_grid, only: &
       GRID_CZ,   &
       GRID_FZ,   &
       GRID_RCDX, &
       GRID_RCDY, &
       GRID_RFDX, &
       GRID_RFDY
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: TOPO_CZX (KA,IA,JA) !< Z coordinate [m] (u,y,interface)
    real(RP) :: TOPO_CZY (KA,IA,JA) !< Z coordinate [m] (x,v,interface)
    real(RP) :: TOPO_CZXY(KA,IA,JA) !< Z coordinate [m] (u,v,interface)
    real(RP) :: TOPO_FZX (KA,IA,JA) !< Z coordinate [m] (u,y,interface)
    real(RP) :: TOPO_FZY (KA,IA,JA) !< Z coordinate [m] (x,v,interface)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! calc Z-coordinate height at staggered position
    do j = 1, JA
    do i = 1, IA-1
       TOPO_CZX(k,i,j) = 0.5D0 * ( TOPO_CZ(k,i+1,j) + TOPO_CZ(k,i,j) )
       TOPO_FZX(k,i,j) = 0.5D0 * ( TOPO_FZ(k,i+1,j) + TOPO_FZ(k,i,j) )
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
       TOPO_CZY(k,i,j) = 0.5D0 * ( TOPO_CZ(k,i,j+1) + TOPO_CZ(k,i,j) )
       TOPO_FZY(k,i,j) = 0.5D0 * ( TOPO_FZ(k,i,j+1) + TOPO_FZ(k,i,j) )
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA-1
       TOPO_CZXY(k,i,j) = 0.25D0 * ( TOPO_CZ(k,i+1,j+1) + TOPO_CZ(k,i+1,j) &
                                   + TOPO_CZ(k,i  ,j+1) + TOPO_CZ(k,i  ,j) )
    enddo
    enddo

    ! G^1/2
    do j = JS, JE
    do i = IS, IE
       ! at (x,y,layer)
       TOPO_GSQRT(1,i,j) = ( TOPO_FZ(1,i,j) - TOPO_CZ(1,i,j) ) &
                         / ( GRID_FZ(1)     - GRID_CZ(1)     )
       do k = 2, KA
          TOPO_GSQRT(k,i,j) = ( TOPO_FZ(k,i,j) - TOPO_FZ(k-1,i,j) ) &
                            / ( GRID_FZ(k)     - GRID_FZ(k-1)     )
       enddo

       ! at (u,y,layer)
       TOPO_GSQRTX(1,i,j) = ( TOPO_FZX(1,i,j) - TOPO_CZX(1,i,j) ) &
                          / ( GRID_FZ (1)     - GRID_CZ (1)     )
       do k = 2, KA
          TOPO_GSQRTX(k,i,j) = ( TOPO_FZX(k,i,j) - TOPO_FZX(k-1,i,j) ) &
                             / ( GRID_FZ (k)     - GRID_FZ (k-1)     )
       enddo

       ! at (x,v,layer)
       TOPO_GSQRTY(1,i,j) = ( TOPO_FZY(1,i,j) - TOPO_CZY(1,i,j) ) &
                          / ( GRID_FZ (1)     - GRID_CZ (1)     )
       do k = 2, KA
          TOPO_GSQRTY(k,i,j) = ( TOPO_FZY(k,i,j) - TOPO_FZY(k-1,i,j) ) &
                             / ( GRID_FZ (k)     - GRID_FZ (k-1)     )
       enddo

       ! at (x,y,interface)
       do k = 1, KA-1
          TOPO_GSQRTZ(k,i,j) = ( TOPO_CZ(k+1,i,j) - TOPO_CZ(k,i,j) ) &
                             / ( GRID_CZ(k+1)     - GRID_CZ(k)     )
       enddo
       TOPO_GSQRTZ(KA,i,j) = ( TOPO_FZ(KA,i,j) - TOPO_CZ(KA,i,j) ) &
                           / ( GRID_FZ(KA)     - GRID_CZ(KA)     )
    enddo
    enddo

    ! Gvector
    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       TOPO_J13(k,i,j) = ( TOPO_CZX(k,i,j) - TOPO_CZX(k,i-1,j  ) ) * GRID_RCDX(i)
       TOPO_J23(k,i,j) = ( TOPO_CZY(k,i,j) - TOPO_CZY(k,i  ,j-1) ) * GRID_RCDY(j)

       TOPO_J13X(k,i,j) = ( TOPO_CZ  (k,i+1,j) - TOPO_CZ  (k,i,j  ) ) * GRID_RFDX(i)
       TOPO_J23X(k,i,j) = ( TOPO_CZXY(k,i  ,j) - TOPO_CZXY(k,i,j-1) ) * GRID_RCDY(j)

       TOPO_J13Y(k,i,j) = ( TOPO_CZXY(k,i,j  ) - TOPO_CZXY(k,i-1,j) ) * GRID_RFDX(i)
       TOPO_J23Y(k,i,j) = ( TOPO_CZ  (k,i,j+1) - TOPO_CZ  (k,i  ,j) ) * GRID_RCDY(j)
    enddo
    enddo
    enddo

    ! fill IHALO & JHALO
    call COMM_vars8( TOPO_GSQRT (:,:,:),  1 )
    call COMM_vars8( TOPO_GSQRTX(:,:,:),  2 )
    call COMM_vars8( TOPO_GSQRTY(:,:,:),  3 )
    call COMM_vars8( TOPO_GSQRTZ(:,:,:),  4 )
    call COMM_vars8( TOPO_J13   (:,:,:),  5 )
    call COMM_vars8( TOPO_J23   (:,:,:),  6 )
    call COMM_vars8( TOPO_J13X  (:,:,:),  7 )
    call COMM_vars8( TOPO_J23X  (:,:,:),  8 )
    call COMM_vars8( TOPO_J13Y  (:,:,:),  9 )
    call COMM_vars8( TOPO_J23Y  (:,:,:), 10 )

    call COMM_wait ( TOPO_GSQRT (:,:,:),  1 )
    call COMM_wait ( TOPO_GSQRTX(:,:,:),  2 )
    call COMM_wait ( TOPO_GSQRTY(:,:,:),  3 )
    call COMM_wait ( TOPO_GSQRTZ(:,:,:),  4 )
    call COMM_wait ( TOPO_J13   (:,:,:),  5 )
    call COMM_wait ( TOPO_J23   (:,:,:),  6 )
    call COMM_wait ( TOPO_J13X  (:,:,:),  7 )
    call COMM_wait ( TOPO_J23X  (:,:,:),  8 )
    call COMM_wait ( TOPO_J13Y  (:,:,:),  9 )
    call COMM_wait ( TOPO_J23Y  (:,:,:), 10 )

    return
  end subroutine TOPO_metrics

end module mod_topography
