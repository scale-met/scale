!-------------------------------------------------------------------------------
!> module atmosphere / grid / cartesC index
!!
!! @par Description
!!          Atmospheric grid Index module for the CartesianC grid
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_grid_cartesC_index
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
  public :: ATMOS_GRID_CARTESC_INDEX_setup
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: ZDIR  = 1
  integer, public, parameter :: XDIR  = 2
  integer, public, parameter :: YDIR  = 3

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

  integer,  public :: I_XYZ = 1 ! at (x,y,z)
  integer,  public :: I_XYW = 2 ! at (x,y,w)
  integer,  public :: I_UYW = 3 ! at (u,y,w)
  integer,  public :: I_XVW = 4 ! at (x,v,w)
  integer,  public :: I_UYZ = 5 ! at (u,y,z)
  integer,  public :: I_XVZ = 6 ! at (x,v,z)
  integer,  public :: I_UVZ = 7 ! at (u,v,z)

  integer,  public :: I_XY  = 1 ! at (x,y)
  integer,  public :: I_UY  = 2 ! at (u,y)
  integer,  public :: I_XV  = 3 ! at (x,v)
  integer,  public :: I_UV  = 4 ! at (u,v)

  integer,  public :: I_FYZ = 1 ! y-z face limiting x-flux
  integer,  public :: I_FXZ = 2 ! x-z face limiting y-flux
  integer,  public :: I_FXY = 3 ! x-y face limiting z-flux

contains

  !-----------------------------------------------------------------------------
  !> setup index
  subroutine ATMOS_GRID_CARTESC_INDEX_setup( &
       KMAX,                &
       IMAXG, JMAXG,        &
       IMAX, JMAX,          &
       KHALO, IHALO, JHALO, &
       IBLOCK, JBLOCK       )
    implicit none
    integer, intent(in), optional :: KMAX
    integer, intent(in), optional :: IMAXG, JMAXG
    integer, intent(in), optional :: IMAX, JMAX
    integer, intent(in), optional :: KHALO, IHALO, JHALO
    integer, intent(in), optional :: IBLOCK, JBLOCK

    call ATMOS_GRID_CARTESC_index_setup_main( &
         KMAX,                &
         IMAXG, JMAXG,        &
         IMAX, JMAX,          &
         KHALO, IHALO, JHALO, &
         IBLOCK, JBLOCK       )

    return
  end subroutine ATMOS_GRID_CARTESC_INDEX_setup

  subroutine ATMOS_GRID_CARTESC_index_setup_main( &
       KMAX_in,                      &
       IMAXG_in, JMAXG_in,           &
       IMAX_in, JMAX_in,             &
       KHALO_in, IHALO_in, JHALO_in, &
       IBLOCK_in, JBLOCK_in          )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    use scale_prc_cartesC, only: &
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
    integer, intent(in), optional :: KMAX_in
    integer, intent(in), optional :: IMAXG_in, JMAXG_in
    integer, intent(in), optional :: IMAX_in, JMAX_in
    integer, intent(in), optional :: KHALO_in, IHALO_in, JHALO_in
    integer, intent(in), optional :: IBLOCK_in, JBLOCK_in

    namelist / PARAM_ATMOS_GRID_CARTESC_INDEX / &
       KMAX,       &
       IMAXG,      &
       JMAXG,      &
       IMAX,       &
       JMAX,       &
       IHALO,      &
       JHALO,      &
       IBLOCK,     &
       JBLOCK

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( present(KMAX_in)  ) KMAX   = KMAX_in
    if ( present(IMAXG_in) ) IMAXG  = IMAXG_in
    if ( present(JMAXG_in) ) JMAXG  = JMAXG_in
    if ( present(IMAX_in)  ) IMAX   = IMAX_in
    if ( present(JMAX_in)  ) JMAX   = JMAX_in
!    if ( present(KHALO_in) ) KHALO  = KHALO_in
    if ( present(IHALO_in) ) IHALO  = IHALO_in
    if ( present(JHALO_in) ) JHALO  = JHALO_in
    if ( present(IBLOCK_in) ) IBLOCK = IBLOCK_in
    if ( present(JBLOCK_in) ) JBLOCK = JBLOCK_in

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[CartesC INDEX] / Categ[ATMOSPHER GRID] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_GRID_CARTESC_INDEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_GRID_CARTESC_index_setup_main",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'Not appropriate names in namelist PARAM_ATMOS_GRID_CARTESC_INDEX. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_GRID_CARTESC_INDEX)



    if ( IMAXG * JMAXG < 0 ) then
       LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'Both IMAXG and JMAXG must set! ', IMAXG, JMAXG
       call PRC_abort
    endif
    if ( IMAX * JMAX < 0 ) then
       LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'Both IMAX and JMAX must set! ', IMAX, JMAX
       call PRC_abort
    endif

    if ( IMAX > 0 .AND. JMAX > 0 ) then
       IMAXG = IMAX * PRC_NUM_X
       JMAXG = JMAX * PRC_NUM_Y
    elseif( IMAXG > 0 .AND. JMAXG > 0 ) then
       IMAX = (IMAXG-1) / PRC_NUM_X + 1
       JMAX = (JMAXG-1) / PRC_NUM_Y + 1

       if ( mod(IMAXG,PRC_NUM_X) > 0 ) then
          LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'number of IMAXG should be divisible by PRC_NUM_X'
          call PRC_abort
!           LOG_INFO("ATMOS_GRID_CARTESC_index_setup_main",*) 'number of IMAXG should be divisible by PRC_NUM_X'
!           LOG_INFO("ATMOS_GRID_CARTESC_index_setup_main",*) 'Small IMAX is used in ranks(X,*)=', PRC_NUM_X-1
!           if ( PRC_2Drank(PRC_myrank,1) == PRC_NUM_X-1 ) then
!              IMAX = IMAXG - IMAX * (PRC_NUM_X-1)
!              LOG_INFO("ATMOS_GRID_CARTESC_index_setup_main",*) 'Small IMAX is used in this rank. IMAX=', IMAX
!           endif
       endif

       if ( mod(JMAXG,PRC_NUM_Y) > 0 ) then
          LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'number of JMAXG should be divisible by PRC_NUM_Y'
          call PRC_abort
!           LOG_INFO("ATMOS_GRID_CARTESC_index_setup_main",*) 'number of JMAXG should be divisible by PRC_NUM_Y'
!           LOG_INFO("ATMOS_GRID_CARTESC_index_setup_main",*) 'Small JMAX is used in ranks(*,Y)=', PRC_NUM_Y-1
!           if ( PRC_2Drank(PRC_myrank,2) == PRC_NUM_Y-1 ) then
!              JMAX = JMAXG - JMAX * (PRC_NUM_Y-1)
!              LOG_INFO("ATMOS_GRID_CARTESC_index_setup_main",*) 'Small JMAX is used in this rank. JMAX=', JMAX
!           endif
       endif
    else
       LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'IMAXG&JMAXG or IMAX&JMAX must set!'
       call PRC_abort
    endif

    if ( IMAX < IHALO ) then
       LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'number of grid size IMAX must >= IHALO! ', IMAX, IHALO
       call PRC_abort
    endif
    if ( JMAX < JHALO ) then
       LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'number of grid size JMAX must >= JHALO! ', JMAX, JHALO
       call PRC_abort
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


    !-- Block size must be divisible
    if    ( mod(IMAX,IBLOCK) > 0 ) then
       LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'number of grid size IMAX must be divisible by IBLOCK! ', IMAX, IBLOCK
       call PRC_abort
    elseif( mod(JMAX,JBLOCK) > 0 ) then
       LOG_ERROR("ATMOS_GRID_CARTESC_index_setup_main",*) 'number of grid size JMAX must be divisible by JBLOCK! ', JMAX, JBLOCK
       call PRC_abort
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

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_index_setup_main",*) 'Atmosphere grid index information '

    ! global
    LOG_INFO_CONT('(1x,3(A,I6))') 'No. of Computational Grid (global)  :', &
                                   KMAX,' x ',IMAXG,' x ',JMAXG
    LOG_INFO_CONT('(1x,2(A,I6))') 'Global index of local grid (X)      :', &
                                   IS_inG," - ",IE_inG
    LOG_INFO_CONT('(1x,2(A,I6))') 'Global index of local grid (Y)      :', &
                                                JS_inG," - ",JE_inG

    ! local
    LOG_NEWLINE
    LOG_INFO_CONT('(1x,3(A,I6))') 'No. of Computational Grid (local)   :', &
                                   KMAX,' x ',IMAX,' x ',JMAX
    LOG_INFO_CONT('(1x,3(A,I6))') 'No. of Grid (including HALO, local) :', &
                                   KA," x ",IA," x ",JA
    LOG_INFO_CONT('(1x,2(A,I6))') 'Local index of inner grid (X)       :', &
                                   ISB," - ",IEB
    LOG_INFO_CONT('(1x,2(A,I6))') 'Local index of inner grid (Y)       :', &
                                                JSB," - ",JEB

    return
  end subroutine ATMOS_GRID_CARTESC_index_setup_main

end module scale_atmos_grid_cartesC_index
