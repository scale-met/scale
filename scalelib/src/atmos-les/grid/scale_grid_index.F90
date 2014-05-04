!-------------------------------------------------------------------------------
!> module grid index
!!
!! @par Description
!!          Grid Index module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-02 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_grid_index
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: KHALO = 2               ! # of halo cells: z
  integer, public, parameter :: IHALO = 2               ! # of halo cells: x
  integer, public, parameter :: JHALO = 2               ! # of halo cells: y

  integer, public, parameter :: ZDIR = 1
  integer, public, parameter :: XDIR = 2
  integer, public, parameter :: YDIR = 3

  integer, public, parameter :: I_DENS = 1
  integer, public, parameter :: I_MOMZ = 2
  integer, public, parameter :: I_MOMX = 3
  integer, public, parameter :: I_MOMY = 4
  integer, public, parameter :: I_RHOT = 5
  integer, public, parameter :: I_QTRC = 6

  integer, public, parameter :: I_BND_VELZ = 1 ! reference velocity (z) [m/s]
  integer, public, parameter :: I_BND_VELX = 2 ! reference velocity (x) [m/s]
  integer, public, parameter :: I_BND_VELY = 3 ! reference velocity (y) [m/s]
  integer, public, parameter :: I_BND_POTT = 4 ! reference potential temperature [K]
  integer, public, parameter :: I_BND_QV   = 5 ! reference water vapor [kg/kg]



#ifdef FIXED_INDEX
  include "inc_index.h"
#else
  integer, public :: KMAX =   -1 ! # of computational cells: z
  integer, public :: IMAX =   -1 ! # of computational cells: x
  integer, public :: JMAX =   -1 ! # of computational cells: y

  integer, public :: IBLOCK = -1 ! block size for cache blocking: x
  integer, public :: JBLOCK = -1 ! block size for cache blocking: y

  integer, public :: KA ! # of z whole cells (local, with HALO)
  integer, public :: IA ! # of x whole cells (local, with HALO)
  integer, public :: JA ! # of y whole cells (local, with HALO)

  integer, public :: KS ! start point of inner domain: z, local
  integer, public :: KE ! end   point of inner domain: z, local
  integer, public :: IS ! start point of inner domain: x, local
  integer, public :: IE ! end   point of inner domain: x, local
  integer, public :: JS ! start point of inner domain: y, local
  integer, public :: JE ! end   point of inner domain: y, local
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  public :: GRID_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine GRID_INDEX_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

#ifndef FIXED_INDEX
    namelist / PARAM_INDEX / &
       KMAX,   &
       IMAX,   &
       JMAX,   &
       IBLOCK, &
       JBLOCK
#endif

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[INDEX]/Categ[COMMON]'

#ifndef FIXED_INDEX
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_INDEX,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_INDEX. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_INDEX)

    KA   = KMAX + KHALO * 2
    IA   = IMAX + IHALO * 2
    JA   = JMAX + JHALO * 2

    KS   = 1    + KHALO
    KE   = KMAX + KHALO
    IS   = 1    + IHALO
    IE   = IMAX + IHALO
    JS   = 1    + JHALO
    JE   = JMAX + JHALO

    IF ( IBLOCK .eq. -1 ) IBLOCK = IMAX
    IF ( JBLOCK .eq. -1 ) JBLOCK = JMAX
#endif

  end subroutine GRID_INDEX_setup

end module scale_grid_index
