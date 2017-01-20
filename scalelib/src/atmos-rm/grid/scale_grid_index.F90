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

  integer, public :: KA          !< # of whole         cells: z, local, with HALO
  integer, public :: IA          !< # of whole         cells: x, local, with HALO
  integer, public :: JA          !< # of whole         cells: y, local, with HALO

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

  ! global size and offset
  integer, public :: IMAXG  = -1 !< # of computational cells: x, global
  integer, public :: JMAXG  = -1 !< # of computational cells: y, global
  integer, public :: IAG         !< # of computational grids
  integer, public :: JAG         !< # of computational grids
  integer, public :: IAGB        !< # of computational grids
  integer, public :: JAGB        !< # of computational grids
  integer, public :: IS_inG      !< start point of the inner domain: cx, global
  integer, public :: IE_inG      !< end   point of the inner domain: cx, global
  integer, public :: JS_inG      !< start point of the inner domain: cy, global
  integer, public :: JE_inG      !< end   point of the inner domain: cy, global
  integer, public :: ISGA        !< start point of the full domain: cx, global
  integer, public :: IEGA        !< end   point of the full domain: cx, global
  integer, public :: JSGA        !< start point of the full domain: cy, global
  integer, public :: JEGA        !< end   point of the full domain: cy, global
  integer, public :: ISGB        !< start point of the inner domain: x, global
  integer, public :: IEGB        !< end   point of the inner domain: x, global
  integer, public :: JSGB        !< start point of the inner domain: y, global
  integer, public :: JEGB        !< end   point of the inner domain: y, global

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
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y, &
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
          write(*,*)                     'xxx number of IMAXG should be divisible by PRC_NUM_X'
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
          write(*,*)                     'xxx number of JMAXG should be divisible by PRC_NUM_Y'
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

    ! array size (global domain)
    IAG   = IMAXG + IHALO * 2
    JAG   = JMAXG + JHALO * 2

    ! horizontal index (global domain)
    IS_inG = IHALO + 1    + PRC_2Drank(PRC_myrank,1) * IMAX
    IE_inG = IHALO + IMAX + PRC_2Drank(PRC_myrank,1) * IMAX
    JS_inG = JHALO + 1    + PRC_2Drank(PRC_myrank,2) * JMAX
    JE_inG = JHALO + JMAX + PRC_2Drank(PRC_myrank,2) * JMAX

    if ( PRC_2Drank(PRC_myrank,1) == 0 ) then
       ISGA = 1
    else
       ISGA = IS_inG
    end if
    if ( PRC_2Drank(PRC_myrank,1) == PRC_NUM_X - 1 ) then
       IEGA = IAG
    else
       IEGA = IE_inG
    end if
    if ( PRC_2Drank(PRC_myrank,2) == 0 ) then
       JSGA = 1
    else
       JSGA = JS_inG
    end if
    if ( PRC_2Drank(PRC_myrank,2) == PRC_NUM_Y - 1 ) then
       JEGA = JAG
    else
       JEGA = JE_inG
    end if

    if ( PRC_PERIODIC_X ) then
       IAGB = IMAXG
       ISGB = IS_inG - IHALO
       IEGB = IE_inG - IHALO
    else
       IAGB = IAG
       if ( PRC_HAS_W ) then
          ISGB = IS_inG
       else ! western boundary
          ISGB = IS_inG - IHALO ! ISGB = 1
       end if
       if ( PRC_HAS_E ) then
          IEGB = IE_inG
       else ! eastern boundary
          IEGB = IE_inG + IHALO
       end if
    end if
    if ( PRC_PERIODIC_Y ) then
       JAGB = JMAXG
       JSGB = JS_inG - JHALO
       JEGB = JE_inG - JHALO
    else
       JAGB = JAG
       if ( PRC_HAS_S ) then
          JSGB = JS_inG
       else ! southern boundary
          JSGB = JS_inG - JHALO ! JSGY = 1
       end if
       if ( PRC_HAS_N ) then
          JEGB = JE_inG
       else ! northern boundary
          JEGB = JE_inG + JHALO
       end if
    end if

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

    ! global
    if( IO_L ) write(IO_FID_LOG,'(1x,3(A,I6))') '*** No. of Computational Grid (global)  :', &
                                                KMAX,' x ',IMAXG,' x ',JMAXG
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))') '*** Global index of local grid (X)      :', &
                                                IS_inG," - ",IE_inG
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))') '*** Global index of local grid (Y)      :', &
                                                JS_inG," - ",JE_inG

    ! local
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,3(A,I6))') '*** No. of Computational Grid (local)   :', &
                                                KMAX,' x ',IMAX,' x ',JMAX
    if( IO_L ) write(IO_FID_LOG,'(1x,3(A,I6))') '*** No. of Grid (including HALO, local) :', &
                                                KA," x ",IA," x ",JA
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))') '*** Local index of inner grid (X)       :', &
                                                ISB," - ",IEB
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))') '*** Local index of inner grid (Y)       :', &
                                                JSB," - ",JEB

    return
  end subroutine GRID_INDEX_setup

end module scale_grid_index
