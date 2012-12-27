!-------------------------------------------------------------------------------
!> module Topography
!!
!! @par Description
!!          Topography module
!!          Terrain-following and vertical metrics
!!
!! @author H.Tomita and SCALE developers
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_topography
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FILECHR, &
     IO_FID_LOG, &
     IO_L
  use gtool_file_h, only: &
     File_HLONG
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TOPO_setup
  public :: TOPO_write
  
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
  real(RP), public, save :: TOPO_Zsfc(1,IA,JA) ! absolute ground height [m]

  real(RP), public, save :: TOPO_CZ(KA,IA,JA) ! Xi center coordinate [m]: z
  real(RP), public, save :: TOPO_FZ(KA,IA,JA) ! Xi face   coordinate [m]: z

  real(RP), public, save :: TOPO_CGSQRT(KA,IA,JA)   ! vertical metrics {G}^1/2
  real(RP), public, save :: TOPO_FGSQRT(KA,IA,JA)   ! vertical metrics {G}^1/2

  real(RP), public, save :: TOPO_Gvec  (KA,IA,JA,2) ! horizontal metrics vector

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
  integer, private, parameter :: XDIR = 1
  integer, private, parameter :: YDIR = 2

  character(len=IO_FILECHR), private, save :: TOPO_IN_BASENAME  = ''
  character(len=IO_FILECHR), private, save :: TOPO_OUT_BASENAME = ''

  character(len=File_HLONG), private, save :: TOPO_OUT_TITLE     = 'SCALE3 TOPOGRAPHY'
  character(len=File_HLONG), private, save :: TOPO_OUT_SOURCE    = 'SCALE-LES ver. 3'
  character(len=File_HLONG), private, save :: TOPO_OUT_INSTITUTE = 'AICS/RIKEN'

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup horizontal&vertical grid
  !-----------------------------------------------------------------------------
  subroutine TOPO_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_TOPO / &
       TOPO_IN_BASENAME,  &
       TOPO_OUT_BASENAME

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

    return
  end subroutine TOPO_setup

  !-----------------------------------------------------------------------------
  !> Read topography
  !-----------------------------------------------------------------------------
  subroutine TOPO_read
    use mod_process, only: &
       PRC_myrank
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use gtool_file, only : &
       FileRead
    implicit none

    real(RP) :: temp(IMAX,JMAX) !> temp file (no HALO)

    character(len=IO_FILECHR) :: bname
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input topography file ***'

    if ( TOPO_IN_BASENAME /= '' ) then

       bname = TOPO_IN_BASENAME

       call FileRead( temp(:,:), bname, 'TOPO', 1, PRC_myrank )

       TOPO_Zsfc(1,IS:IE,JS:JE) = temp(1:IMAX,1:JMAX)

       ! fill IHALO & JHALO
       call COMM_vars8( TOPO_Zsfc(1,:,:), 1 )
       call COMM_wait ( TOPO_Zsfc(1,:,:), 1 )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** topography file is not specified.'

       TOPO_Zsfc(:,:,:) = 0.D0
    endif

    return
  end subroutine TOPO_read

  !-----------------------------------------------------------------------------
  !> Write topography
  !-----------------------------------------------------------------------------
  subroutine TOPO_write
    use mod_process, only: &
       PRC_master, &
       PRC_myrank
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate, &
       FileAddVariable, &
       FilePutAxis, &
       FileWrite, &
       FileClose
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_grid, only :  &
       GRID_CX, &
       GRID_CY
    implicit none

    real(RP) :: temp(IMAX,JMAX) !> temp file (no HALO)

    character(len=IO_FILECHR) :: bname

    integer :: dtype
    integer :: fid, vid
    !---------------------------------------------------------------------------

    if ( TOPO_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output topography file ***'

       ! fill IHALO & JHALO
       call COMM_vars8( TOPO_Zsfc(1,:,:), 1 )
       call COMM_wait ( TOPO_Zsfc(1,:,:), 1 )

       temp(1:IMAX,1:JMAX) = TOPO_Zsfc(1,IS:IE,JS:JE)

       bname = TOPO_OUT_BASENAME

       call FileCreate( fid,                                     & ! [OUT]
                        bname,                                   & ! [IN]
                        TOPO_OUT_TITLE,                          & ! [IN]
                        TOPO_OUT_SOURCE,                         & ! [IN]
                        TOPO_OUT_INSTITUTE,                      & ! [IN]
                        (/'x','y'/), (/IMAX,JMAX/), (/'X','Y'/), & ! [IN]
                        (/'m','m'/), (/File_REAL4,File_REAL4/),  & ! [IN]
                        PRC_master, PRC_myrank                   ) ! [IN]

       call FilePutAxis( fid, 'x', GRID_CX(IS:IE) )
       call FilePutAxis( fid, 'y', GRID_CY(JS:JE) )

       if    ( RP == 8 ) then
          dtype = File_REAL8
       elseif( RP == 4 ) then
          dtype = File_REAL4
       endif

       call FileAddVariable( vid,                              & ! [OUT]
                             fid, 'TOPO', 'Topography', '[m]', & ! [IN]
                             (/'x','y'/), dtype                ) ! [IN]

       call FileWrite( vid, temp(:,:), NOWSEC, NOWSEC )

       call FileClose( fid )

    endif

    return
  end subroutine TOPO_write

  !-----------------------------------------------------------------------------
  !> Convert Z to Xi coordinate
  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------
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
