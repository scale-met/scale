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
  integer, public, parameter :: ZDIR  = 1
  integer, public, parameter :: XDIR  = 2
  integer, public, parameter :: YDIR  = 3

#ifdef _FIXEDINDEX_
  include "inc_index.h"
  include "inc_index_common.h"
#else
  integer, public :: KMAX   = -1 !< # of computational cells: z
  integer, public :: IMAX   = -1 !< # of computational cells: x
  integer, public :: JMAX   = -1 !< # of computational cells: y

  integer, public :: IBLOCK = -1 !< block size for cache blocking: x
  integer, public :: JBLOCK = -1 !< block size for cache blocking: y

  integer, public, parameter :: KHALO  = 2  !< # of halo cells: z
  integer, public            :: IHALO  = 2  !< # of halo cells: x
  integer, public            :: JHALO  = 2  !< # of halo cells: y

  integer, public :: KA          !< # of z whole cells (local, with HALO)
  integer, public :: IA          !< # of x whole cells (local, with HALO)
  integer, public :: JA          !< # of y whole cells (local, with HALO)

  integer, public :: KS          !< start point of inner domain: z, local
  integer, public :: KE          !< end   point of inner domain: z, local
  integer, public :: IS          !< start point of inner domain: x, local
  integer, public :: IE          !< end   point of inner domain: x, local
  integer, public :: JS          !< start point of inner domain: y, local
  integer, public :: JE          !< end   point of inner domain: y, local

  integer, public :: KIJMAX = -1 !< # of computational cells: z*x*y
#endif

  ! indices considering boundary
  integer, public :: IMAXB
  integer, public :: JMAXB
  integer, public :: ISB
  integer, public :: IEB
  integer, public :: JSB
  integer, public :: JEB
  integer, public :: IEH         !< end   point of inner domain: x, local (half level)
  integer, public :: JEH         !< end   point of inner domain: y, local (half level)

  integer, public :: ISG         !< start point of the inner domain: x, global
  integer, public :: IEG         !< end   point of the inner domain: x, global
  integer, public :: JSG         !< start point of the inner domain: y, global
  integer, public :: JEG         !< end   point of the inner domain: y, global

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
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank,  &
       PRC_NUM_X,   &
       PRC_NUM_Y,   &
       PRC_HAS_W,   &
       PRC_HAS_E,   &
       PRC_HAS_S,   &
       PRC_HAS_N
    implicit none

#ifndef _FIXEDINDEX_
    namelist / PARAM_INDEX / &
       KMAX,   &
       IMAX,   &
       JMAX,   &
       IHALO,  &
       JHALO,  &
       IBLOCK, &
       JBLOCK
#endif

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRID_INDEX] / Categ[ATMOS-RM GRID] / Origin[SCALElib]'

#ifdef _FIXEDINDEX_
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

    if ( IMAX < IHALO ) then
       write(*,*) 'xxx number of grid size IMAX must >= IHALO! ', IMAX, IHALO
       call PRC_MPIstop
    endif
    if ( JMAX < JHALO ) then
       write(*,*) 'xxx number of grid size JMAX must >= JHALO! ', JMAX, JHALO
       call PRC_MPIstop
    endif

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

    KIJMAX = KMAX * IMAX * JMAX
#endif

    !-- Block size must be divisible
    if    ( mod(IMAX,IBLOCK) > 0 ) then
       write(*,*) 'xxx number of grid size IMAX must be divisible by IBLOCK! ', IMAX, IBLOCK
       call PRC_MPIstop
    elseif( mod(JMAX,JBLOCK) > 0 ) then
       write(*,*) 'xxx number of grid size JMAX must be divisible by JBLOCK! ', JMAX, JBLOCK
       call PRC_MPIstop
    endif

    ! horizontal index (global domain)
    ISG = IHALO + 1    + PRC_2Drank(PRC_myrank,1) * IMAX
    IEG = IHALO + IMAX + PRC_2Drank(PRC_myrank,1) * IMAX
    JSG = JHALO + 1    + PRC_2Drank(PRC_myrank,2) * JMAX
    JEG = JHALO + JMAX + PRC_2Drank(PRC_myrank,2) * JMAX

    ! index considering boundary region
    IMAXB = IMAX
    JMAXB = JMAX
    ISB = IS
    IEB = IE
    JSB = JS
    JEB = JE
    IEH = IE
    JEH = JE
    if ( .NOT. PRC_HAS_W ) then
       IMAXB = IMAXB + IHALO
       ISB = 1
    endif
    if ( .NOT. PRC_HAS_E ) then
       IMAXB = IMAXB + IHALO
       IEB = IA
       IEH = IE - 1
    endif
    if ( .NOT. PRC_HAS_S ) then
       JMAXB = JMAXB + JHALO
       JSB = 1
    endif
    if ( .NOT. PRC_HAS_N ) then
       JMAXB = JMAXB + JHALO
       JEB = JA
       JEH = JE - 1
    endif

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
