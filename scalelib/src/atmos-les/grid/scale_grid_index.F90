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
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GRID_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: KHALO = 2 ! # of halo cells: z
  integer, public, parameter :: IHALO = 2 ! # of halo cells: x
  integer, public, parameter :: JHALO = 2 ! # of halo cells: y

  integer, public, parameter :: ZDIR  = 1
  integer, public, parameter :: XDIR  = 2
  integer, public, parameter :: YDIR  = 3

#ifdef FIXED_INDEX
  include "inc_index.h"
#else
  integer, public            :: KMAX =   -1 !< # of computational cells: z
  integer, public            :: IMAX =   -1 !< # of computational cells: x
  integer, public            :: JMAX =   -1 !< # of computational cells: y

  integer, public            :: IBLOCK = -1 !< block size for cache blocking: x
  integer, public            :: JBLOCK = -1 !< block size for cache blocking: y

  integer, public            :: KA          !< # of z whole cells (local, with HALO)
  integer, public            :: IA          !< # of x whole cells (local, with HALO)
  integer, public            :: JA          !< # of y whole cells (local, with HALO)

  integer, public            :: KS          !< start point of inner domain: z, local
  integer, public            :: KE          !< end   point of inner domain: z, local
  integer, public            :: IS          !< start point of inner domain: x, local
  integer, public            :: IE          !< end   point of inner domain: x, local
  integer, public            :: JS          !< start point of inner domain: y, local
  integer, public            :: JE          !< end   point of inner domain: y, local
#endif

  integer, public            :: ISG         !< start point of the inner domain: x, global
  integer, public            :: IEG         !< end   point of the inner domain: x, global
  integer, public            :: JSG         !< start point of the inner domain: y, global
  integer, public            :: JEG         !< end   point of the inner domain: y, global

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
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
       PRC_MPIstop, &
       PRC_myrank,  &
       PRC_2Drank,  &
       PRC_NUM_X,   &
       PRC_NUM_Y
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
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRID_INDEX] / Categ[ATMOS-LES GRID] / Origin[SCALElib]'

#ifdef FIXED_INDEX
    if( IO_L ) write(IO_FID_LOG,*) '*** No namelists.'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** fixed index mode'
#else
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_INDEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_INDEX. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_INDEX)

    KA   = KMAX + KHALO * 2
    IA   = IMAX + IHALO * 2
    JA   = JMAX + JHALO * 2

    KS   = 1    + KHALO
    KE   = KMAX + KHALO
    IS   = 1    + IHALO
    IE   = IMAX + IHALO
    JS   = 1    + JHALO
    JE   = JMAX + JHALO

    if( IBLOCK == -1 ) IBLOCK = IMAX
    if( JBLOCK == -1 ) JBLOCK = JMAX
#endif

    ! horizontal index (global domain)
    ISG = IHALO + 1    + PRC_2Drank(PRC_myrank,1) * IMAX
    IEG = IHALO + IMAX + PRC_2Drank(PRC_myrank,1) * IMAX
    JSG = JHALO + 1    + PRC_2Drank(PRC_myrank,2) * JMAX
    JEG = JHALO + JMAX + PRC_2Drank(PRC_myrank,2) * JMAX

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Atmosphere grid index information ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)') '*** No. of Computational Grid (global)  :', &
                                                       KMAX,           ' x ',                       &
                                                       IMAX*PRC_NUM_X, ' x ',                       &
                                                       JMAX*PRC_NUM_Y
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)') '*** No. of Computational Grid (local)   :', &
                                                       KMAX,' x ',IMAX,' x ',JMAX
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)') '*** No. of Grid (including HALO, local) :', &
                                                       KA," x ",IA," x ",JA
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6)')      '*** Global index of local grid (X)      :', &
                                                       ISG," - ",IEG
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6)')      '*** Global index of local grid (Y)      :', &
                                                       JSG," - ",JEG

  end subroutine GRID_INDEX_setup

end module scale_grid_index
