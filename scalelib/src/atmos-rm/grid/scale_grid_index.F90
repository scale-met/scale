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
  integer, public :: KMAX   = -1 !< # of computational cells: z, local
  integer, public :: IMAX   = -1 !< # of computational cells: x, local
  integer, public :: JMAX   = -1 !< # of computational cells: y, local

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

  integer, public :: IMAXG  = -1 !< # of computational cells: x, global
  integer, public :: JMAXG  = -1 !< # of computational cells: y, global

  integer, public :: ISG         !< start point of the inner domain: x, global
  integer, public :: IEG         !< end   point of the inner domain: x, global
  integer, public :: JSG         !< start point of the inner domain: y, global
  integer, public :: JEG         !< end   point of the inner domain: y, global

  ! indices considering boundary
  integer, public :: IMAXB
  integer, public :: JMAXB
  integer, public :: ISB
  integer, public :: IEB
  integer, public :: JSB
  integer, public :: JEB
  integer, public :: IEH         !< end   point of inner domain: x, local (half level)
  integer, public :: JEH         !< end   point of inner domain: y, local (half level)

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
       IMAXG,  &
       JMAXG,  &
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

    IMAXG = IMAX * PRC_NUM_Y
    JMAXG = JMAX * PRC_NUM_Y
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

    if ( IMAXG * JMAXG < 0 ) then
       write(*,*) 'xxx Both IMAXG and JMAXG must set! ', IMAXG, JMAXG
       call PRC_MPIstop
    endif
    if ( IMAX * JMAX < 0 ) then
       write(*,*) 'xxx Both IMAX and JMAX must set! ', IMAX, JMAX
       call PRC_MPIstop
    endif

    if ( IMAX > 0 .AND. JMAX > 0 ) then
       IMAXG = IMAX * PRC_NUM_X
       JMAXG = JMAX * PRC_NUM_Y
    elseif( IMAXG > 0 .AND. JMAXG > 0 ) then
       IMAX = (IMAXG-1) / PRC_NUM_X + 1
       JMAX = (JMAXG-1) / PRC_NUM_Y + 1

       if ( mod(IMAXG,PRC_NUM_X) > 0 ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx number of IMAXG should be divisible by PRC_NUM_X'
          if( IO_L ) write(*         ,*) 'xxx number of IMAXG should be divisible by PRC_NUM_X'
          call PRC_MPIstop
!           if( IO_L ) write(IO_FID_LOG,*) '*** number of IMAXG should be divisible by PRC_NUM_X'
!           if( IO_L ) write(IO_FID_LOG,*) '*** Small IMAX is used in ranks(X,*)=', PRC_NUM_X-1
!           if ( PRC_2Drank(PRC_myrank,1) == PRC_NUM_X-1 ) then
!              IMAX = IMAXG - IMAX * (PRC_NUM_X-1)
!              if( IO_L ) write(IO_FID_LOG,*) '*** Small IMAX is used in this rank. IMAX=', IMAX
!           endif
       endif

       if ( mod(JMAXG,PRC_NUM_Y) > 0 ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx number of JMAXG should be divisible by PRC_NUM_Y'
          if( IO_L ) write(*         ,*) 'xxx number of JMAXG should be divisible by PRC_NUM_Y'
          call PRC_MPIstop
!           if( IO_L ) write(IO_FID_LOG,*) '*** number of JMAXG should be divisible by PRC_NUM_Y'
!           if( IO_L ) write(IO_FID_LOG,*) '*** Small JMAX is used in ranks(*,Y)=', PRC_NUM_Y-1
!           if ( PRC_2Drank(PRC_myrank,2) == PRC_NUM_Y-1 ) then
!              JMAX = JMAXG - JMAX * (PRC_NUM_Y-1)
!              if( IO_L ) write(IO_FID_LOG,*) '*** Small JMAX is used in this rank. JMAX=', JMAX
!           endif
       endif
    else
       write(*,*) 'xxx IMAXG&JMAXG or IMAX&JMAX must set!'
       call PRC_MPIstop
    endif

    if ( IMAX < IHALO ) then
       write(*,*) 'xxx number of grid size IMAX must >= IHALO! ', IMAX, IHALO
       call PRC_MPIstop
    endif
    if ( JMAX < JHALO ) then
       write(*,*) 'xxx number of grid size JMAX must >= JHALO! ', JMAX, JHALO
       call PRC_MPIstop
    endif

    KA = KMAX + KHALO * 2
    IA = IMAX + IHALO * 2
    JA = JMAX + JHALO * 2

    KS = 1    + KHALO
    KE = KMAX + KHALO
    IS = 1    + IHALO
    IE = IMAX + IHALO
    JS = 1    + JHALO
    JE = JMAX + JHALO

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
    ISB   = IS
    IEB   = IE
    JSB   = JS
    JEB   = JE
    IEH   = IE
    JEH   = JE

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
    if( IO_L ) write(IO_FID_LOG,'(1x,3(A,I6))') '*** No. of Computational Grid (global)  :', &
                                                KMAX,' x ',IMAXG,' x ',JMAXG

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,3(A,I6))') '*** No. of Computational Grid (local)   :', &
                                                KMAX,' x ',IMAX,' x ',JMAX
    if( IO_L ) write(IO_FID_LOG,'(1x,3(A,I6))') '*** No. of Grid (including HALO, local) :', &
                                                KA," x ",IA," x ",JA
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))') '*** Local index of inner grid (X)       :', &
                                                ISB," - ",IEB
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))') '*** Local index of inner grid (Y)       :', &
                                                JSB," - ",JEB

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))') '*** Global index of local grid (X)      :', &
                                                ISG," - ",IEG
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))') '*** Global index of local grid (Y)      :', &
                                                JSG," - ",JEG

  end subroutine GRID_INDEX_setup

end module scale_grid_index
