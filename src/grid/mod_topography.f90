!-------------------------------------------------------------------------------
!> module Topography
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

  real(RP), public, save :: TOPO_CZ    (KA,IA,JA)   !< Xi center coordinate [m]: z
  real(RP), public, save :: TOPO_FZ    (KA,IA,JA)   !< Xi face   coordinate [m]: z

  real(RP), public, save :: TOPO_CGSQRT(KA,IA,JA)   !< vertical metrics {G}^1/2
  real(RP), public, save :: TOPO_FGSQRT(KA,IA,JA)   !< vertical metrics {G}^1/2

  real(RP), public, save :: TOPO_Gvec  (KA,IA,JA,2) !< horizontal metrics vector

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TOPO_read
  private :: TOPO_Xi
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
    call TOPO_Xi
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

       call FILEIO_write( TOPO_Zsfc(1,:,:),  TOPO_OUT_BASENAME, TOPO_OUT_TITLE, & ! [IN]
                          'TOPO', 'Topography', 'm', 'XY',      TOPO_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine TOPO_write

  !-----------------------------------------------------------------------------
  !> Convert Z to Xi coordinate
  subroutine TOPO_Xi
    use mod_grid, only: &
       GRID_CZ, &
       GRID_FZ
    implicit none

    real(RP) :: Htop

    integer :: k, i, j
    !---------------------------------------------------------------------------

    Htop = GRID_FZ(KE) - GRID_FZ(KS-1)

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       TOPO_CZ(k,i,j) = Htop * ( GRID_CZ(k) - TOPO_Zsfc(1,i,j) ) &
                             / ( Htop       - TOPO_Zsfc(1,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       TOPO_FZ(k,i,j) = Htop * ( GRID_FZ(k) - TOPO_Zsfc(1,i,j) ) &
                             / ( Htop       - TOPO_Zsfc(1,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine TOPO_Xi

  !-----------------------------------------------------------------------------
  !> Calculate G^1/2 & Gvector
  subroutine TOPO_metrics
    use mod_grid, only: &
       GRID_CZ,   &
       GRID_FZ,   &
       GRID_RFDX, &
       GRID_RFDY
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), parameter :: FACT_N =   9.0_RP / 8.0_RP
    real(RP), parameter :: FACT_F = - 1.0_RP / 8.0_RP / 3.0_RP

    real(RP) :: grad

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! G^1/2
    do j = 1, JA
    do i = 1, IA
       TOPO_CGSQRT(1,i,j) = ( TOPO_FZ(1,i,j) - TOPO_CZ(1,i,j) ) &
                          / ( GRID_FZ(1)     - GRID_CZ(1)     )
       do k = 2, KA
          TOPO_CGSQRT(k,i,j) = ( TOPO_FZ(k,i,j) - TOPO_FZ(k-1,i,j) ) &
                             / ( GRID_FZ(k)     - GRID_FZ(k-1)     )
       enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       do k = 1, KA-1
          TOPO_FGSQRT(k,i,j) = ( TOPO_CZ(k+1,i,j) - TOPO_CZ(k,i,j) ) &
                             / ( GRID_CZ(k+1)     - GRID_CZ(k)     )
       enddo
       TOPO_FGSQRT(KA,i,j) = ( TOPO_FZ(KA,i,j) - TOPO_CZ(KA,i,j) ) &
                           / ( GRID_FZ(KA)     - GRID_CZ(KA)     )
    enddo
    enddo

    ! Gvector
    do j = 1,  JA
    do i = IS, IE
    do k = 1,  KA
       grad = ( FACT_N * ( TOPO_CZ(k,i+1,j)-TOPO_CZ(k,i  ,j) ) &
              + FACT_F * ( TOPO_CZ(k,i+2,j)-TOPO_CZ(k,i-1,j) ) ) * GRID_RFDX(i)

       TOPO_Gvec(k,i,j,XDIR) = - grad / TOPO_FGSQRT(k,i,j)
    enddo
    enddo
    enddo


    do j = JS, JE
    do i = 1,  IA
    do k = 1,  KA
       grad = ( FACT_N * ( TOPO_CZ(k,i,j+1)-TOPO_CZ(k,i,j  ) ) &
              + FACT_F * ( TOPO_CZ(k,i,j+2)-TOPO_CZ(k,i,j-1) ) ) * GRID_RFDY(j)

       TOPO_Gvec(k,i,j,YDIR) = - grad / TOPO_FGSQRT(k,i,j)
    enddo
    enddo
    enddo

    ! fill IHALO & JHALO
    call COMM_vars8( TOPO_Gvec(:,:,:,XDIR), 1 )
    call COMM_vars8( TOPO_Gvec(:,:,:,YDIR), 2 )
    call COMM_wait ( TOPO_Gvec(:,:,:,XDIR), 1 )
    call COMM_wait ( TOPO_Gvec(:,:,:,YDIR), 2 )

    return
  end subroutine TOPO_metrics

end module mod_topography
