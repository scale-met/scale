!-------------------------------------------------------------------------------
!> module HISTORY
!!
!! @par Description
!!          History output module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-05 (H.Yashiro)   [new]
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-06-11 (S.Nishizawa) [mod] use gtool_history
!!
!<
!-------------------------------------------------------------------------------
module scale_history
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_land_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: HIST_setup
  public :: HIST_switch
  public :: HIST_setpres
  public :: HIST_reg
  public :: HIST_query
  public :: HIST_put
  public :: HIST_in
  public :: HIST_get
  public :: HIST_write

  interface HIST_put
     module procedure HIST_put_0D
     module procedure HIST_put_1D
     module procedure HIST_put_2D
     module procedure HIST_put_3D
  end interface HIST_put

  interface HIST_in
     module procedure HIST_in_0D
     module procedure HIST_in_1D
     module procedure HIST_in_2D
     module procedure HIST_in_3D
  end interface HIST_in

  interface HIST_get
     module procedure HIST_get_1D
     module procedure HIST_get_2D
     module procedure HIST_get_3D
  end interface HIST_get

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: HIST_put_axes

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter   :: I_MODEL = 0            !< model coordinate
  integer,                private, parameter   :: I_Z     = 1            !< z coordinate
  integer,                private, parameter   :: I_PRES  = 2            !< pressure coordinate

  integer,                private              :: HIST_item_limit
  integer,                private              :: HIST_item_count        !< number of the history item
  character(len=H_SHORT), private, allocatable :: HIST_item   (:)        !< name   of the history item
  integer,                private, allocatable :: HIST_variant(:)        !< number of the variants      for each history item
  integer,                private, allocatable :: HIST_zcoord (:,:)      !< vertical interpolation type for each variant of history item

  integer,                private, parameter   :: HIST_PRES_nlim   = 300 !< limit  of index size for pressure layer
  integer,                private              :: HIST_PRES_nlayer =  -1 !< Number of pressure layer
  real(RP),               private, allocatable :: HIST_PRES_val(:)       !< pressure level to output [hPa]

  logical,                private              :: enabled
  integer,                private              :: im,  jm,  km
  integer,                private              :: ims, ime
  integer,                private              :: jms, jme

  logical, private :: HIST_BND = .false.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine HIST_setup
    use gtool_history, only: &
       HistoryInit
    use scale_process, only: &
       PRC_MPIstop,    &
       PRC_masterrank, &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    use scale_time, only: &
       TIME_NOWDATE,     &
       TIME_NOWMS,       &
       TIME_DTSEC,       &
       TIME_STARTDAYSEC, &
       TIME_OFFSET_YEAR
    use scale_interpolation, only: &
       INTERP_setup_pres
    implicit none

    character(len=H_MID) :: HISTORY_H_TITLE = 'SCALE-RM HISTORY OUTPUT' !< title of the output file
    character(len=H_MID) :: HISTORY_T_SINCE
    real(RP)             :: HIST_PRES(HIST_PRES_nlim)                   !< pressure level to output [hPa]

    NAMELIST / PARAM_HIST / &
       HIST_PRES_nlayer, &
       HIST_PRES,        &
       HIST_BND

    integer  :: HIST_variant_limit

    real(DP) :: start_daysec
    integer  :: rankidx(2)
    integer  :: ierr
    integer  :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[HISTORY] / Categ[ATMOS-RM IO] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_HIST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_HIST. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_HIST)

    ! check pressure coordinate
    if ( HIST_PRES_nlayer > 0 ) then
       if ( HIST_PRES_nlayer > 100 ) then
          write(*,*) 'xxx number of layers of pressure is larger the KMAX'
          call PRC_MPIstop
       end if

       do k = 1, HIST_PRES_nlayer
          if ( HIST_PRES(k) <= 0.0_RP ) then
             write(*,*) 'xxx Invalid value found in pressure coordinate! (k,value)=', k, HIST_PRES(k)
             call PRC_MPIstop
          elseif ( HIST_PRES(k+1) >= HIST_PRES(k) ) then
             write(*,*) 'xxx The value of pressure coordinate must be descending order! ', &
             '(k,value[k],value[k+1])=', k, HIST_PRES(k), HIST_PRES(k+1)
             call PRC_MPIstop
          endif
       enddo
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** HIST_PRES_nlayer is not set.'
       if( IO_L ) write(IO_FID_LOG,*) '*** Output with pressure coordinate is disabled'
    endif

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start_daysec = TIME_STARTDAYSEC
    if ( TIME_NOWDATE(1) > 0 ) then
       write(HISTORY_T_SINCE,'(I4.4,5(A1,I2.2))')      TIME_NOWDATE(1), &
                                                  '-', TIME_NOWDATE(2), &
                                                  '-', TIME_NOWDATE(3), &
                                                  ' ', TIME_NOWDATE(4), &
                                                  ':', TIME_NOWDATE(5), &
                                                  ':', TIME_NOWDATE(6)
       start_daysec = TIME_NOWMS
    else
       HISTORY_T_SINCE = ''
    endif

    if ( HIST_BND ) then
       im  = IMAXB
       jm  = JMAXB
       ims = ISB
       ime = IEB
       jms = JSB
       jme = JEB
    else
       im  = IMAX
       jm  = JMAX
       ims = IS
       ime = IE
       jms = JS
       jme = JE
    end if

    km = max( LKMAX, UKMAX, KMAX )

    call HistoryInit( HIST_item_limit,                  & ! [OUT]
                      HIST_variant_limit,               & ! [OUT]
                      im, jm, km,                       & ! [IN]
                      PRC_masterrank,                   & ! [IN]
                      PRC_myrank,                       & ! [IN]
                      rankidx(:),                       & ! [IN]
                      HISTORY_H_TITLE,                  & ! [IN]
                      H_SOURCE,                         & ! [IN]
                      H_INSTITUTE,                      & ! [IN]
                      time_start     = start_daysec,    & ! [IN]
                      time_interval  = TIME_DTSEC,      & ! [IN]
                      time_since     = HISTORY_T_SINCE, & ! [IN]
                      default_zcoord = 'model',         & ! [IN]
                      namelist_fid   = IO_FID_CONF      ) ! [IN]

    HIST_item_count = 0
    if ( HIST_item_limit > 0 ) then
       allocate( HIST_item   (HIST_item_limit)                    )
       allocate( HIST_variant(HIST_item_limit)                    )
       allocate( HIST_zcoord (HIST_item_limit,HIST_variant_limit) )
       HIST_item   (:)   = ''
       HIST_variant(:)   = 0
       HIST_zcoord (:,:) = 0
    endif

    if ( HIST_PRES_nlayer > 0 ) then
       allocate( HIST_PRES_val(HIST_PRES_nlayer) )

       do k = 1, HIST_PRES_nlayer
          HIST_PRES_val(k) = HIST_PRES(k) * 100.0_RP ! [hPa->Pa]
       enddo

       call INTERP_setup_pres( HIST_PRES_nlayer ) ! [IN]
    endif

    call HIST_put_axes

    enabled = .true.

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_setup

  !-----------------------------------------------------------------------------
  !> Put axis coordinate to history file
  !  only register the axis and coordinate variables into internal buffers
  !  The actual write happens later when calling HIST_write
  subroutine HIST_put_axes
    use gtool_history, only: &
       HistoryPutAxis,  &
       HistoryPutAssociatedCoordinates
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_grid, only: &
       GRID_CZ,    &
       GRID_CX,    &
       GRID_CY,    &
       GRID_FZ,    &
       GRID_FX,    &
       GRID_FY,    &
       GRID_CDZ,   &
       GRID_CDX,   &
       GRID_CDY,   &
       GRID_FDZ,   &
       GRID_FDX,   &
       GRID_FDY,   &
       GRID_CBFZ,  &
       GRID_CBFX,  &
       GRID_CBFY,  &
       GRID_FBFZ,  &
       GRID_FBFX,  &
       GRID_FBFY,  &
       GRID_CXG,   &
       GRID_CYG,   &
       GRID_FXG,   &
       GRID_FYG,   &
       GRID_CDXG,  &
       GRID_CDYG,  &
       GRID_FDXG,  &
       GRID_FDYG,  &
       GRID_CBFXG, &
       GRID_CBFYG, &
       GRID_FBFXG, &
       GRID_FBFYG
    use scale_land_grid, only: &
       GRID_LCZ, &
       GRID_LFZ, &
       GRID_LCDZ
    use scale_urban_grid, only: &
       GRID_UCZ, &
       GRID_UFZ, &
       GRID_UCDZ
    use scale_grid_real, only: &
       REAL_CZ,   &
       REAL_FZ,   &
       REAL_LON,  &
       REAL_LONX, &
       REAL_LONY, &
       REAL_LONXY, &
       REAL_LAT,  &
       REAL_LATX, &
       REAL_LATY, &
       REAL_LATXY
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_landuse, only: &
       LANDUSE_frac_land
    use scale_process, only: &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none

    real(RP)         :: AXIS     (im,jm,KMAX)
    character(len=2) :: AXIS_name(3)

    integer :: k, i, j
    integer :: rankidx(2)
    integer :: start(3)
    integer :: startX, startY, startZ
    integer :: XAG, YAG
    !---------------------------------------------------------------------------

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    ! For parallel I/O, some variables are written by a subset of processes.
    ! 1. Only PRC_myrank 0 writes all z axes
    ! 2. Only south-most processes (rankidx(2) == 0) write x axes
    !         rankidx(1) == 0           writes west HALO
    !         rankidx(1) == PRC_NUM_X-1 writes east HALO
    !         others                    writes without HALO
    ! 3. Only west-most processes (rankidx(1) == 0) write y axes
    !         rankidx(1) == 0           writes south HALO
    !         rankidx(1) == PRC_NUM_Y-1 writes north HALO
    !         others                    writes without HALO

    if ( HIST_BND ) then
       startX = ISGB   ! global subarray starting index
       startY = JSGB   ! global subarray starting index
       XAG    = IAGB
       YAG    = JAGB
    else
       startX = 1 + PRC_2Drank(PRC_myrank,1) * IMAX    ! no IHALO
       startY = 1 + PRC_2Drank(PRC_myrank,2) * JMAX    ! no JHALO
       XAG    = IMAXG
       YAG    = JMAXG
    end if
    startZ = 1
    if ( rankidx(2) .GT. 0 ) startX = -1   ! only south-most processes write
    if ( rankidx(1) .GT. 0 ) startY = -1   ! only west-most processes write
    if ( PRC_myrank .GT. 0 ) startZ = -1   ! only rank 0 writes Z axes

    ! for the shared-file I/O method, the axes are global (gsize)
    ! for one-file-per-process I/O method, the axes size is equal to the local buffer size
    call HistoryPutAxis( 'x',   'X',               'm', 'x',   GRID_CX(ims:ime), gsize=XAG,  start=startX )
    call HistoryPutAxis( 'y',   'Y',               'm', 'y',   GRID_CY(jms:jme), gsize=YAG,  start=startY )
    call HistoryPutAxis( 'z',   'Z',               'm', 'z',   GRID_CZ(KS:KE),   gsize=KMAX, start=startZ )
    call HistoryPutAxis( 'xh',  'X (half level)',  'm', 'xh',  GRID_FX(ims:ime), gsize=XAG,  start=startX )
    call HistoryPutAxis( 'yh',  'Y (half level)',  'm', 'yh',  GRID_FY(jms:jme), gsize=YAG,  start=startY )
    call HistoryPutAxis( 'zh',  'Z (half level)',  'm', 'zh',  GRID_FZ(KS:KE),   gsize=KMAX, start=startZ )

    if ( HIST_PRES_nlayer > 0 ) then
       call HistoryPutAxis( 'pressure', 'Pressure', 'hPa', 'pressure',         &
                            HIST_PRES_val(:)/100.0_RP, down=.true. )
    endif

    call HistoryPutAxis( 'lz',  'LZ',              'm', 'lz',  GRID_LCZ(LKS:LKE), down=.true., gsize=LKMAX, start=startZ )
    call HistoryPutAxis( 'lzh', 'LZ (half level)', 'm', 'lzh', GRID_LFZ(LKS:LKE), down=.true., gsize=LKMAX, start=startZ )

    call HistoryPutAxis( 'uz',  'UZ',              'm', 'uz',  GRID_UCZ(UKS:UKE), down=.true., gsize=UKMAX, start=startZ )
    call HistoryPutAxis( 'uzh', 'UZ (half level)', 'm', 'uzh', GRID_UFZ(UKS:UKE), down=.true., gsize=UKMAX, start=startZ )

    ! axes below always include halos when written to file regardless of PRC_PERIODIC_X/PRC_PERIODIC_Y
    call HistoryPutAxis( 'CZ',  'Atmos Grid Center Position Z', 'm', 'CZ',  GRID_CZ,  gsize=KA,    start=startZ )
    if ( IO_PNETCDF ) then
       call HistoryPutAxis( 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  GRID_CXG, gsize=IAG,   start=startZ )
       call HistoryPutAxis( 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  GRID_CYG, gsize=JAG,   start=startZ )
    else
       call HistoryPutAxis( 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  GRID_CX )
       call HistoryPutAxis( 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  GRID_CY )
    end if
    call HistoryPutAxis( 'FZ',  'Atmos Grid Face Position Z',   'm', 'FZ',  GRID_FZ,  gsize=KA+1,  start=startZ )
    if ( IO_PNETCDF ) then
       call HistoryPutAxis( 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  GRID_FXG, gsize=IAG+1, start=startZ )
       call HistoryPutAxis( 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  GRID_FYG, gsize=JAG+1, start=startZ )
    else
       call HistoryPutAxis( 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  GRID_FX )
       call HistoryPutAxis( 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  GRID_FY )
    end if

    call HistoryPutAxis( 'CDZ',  'Grid Cell length Z', 'm', 'CZ',  GRID_CDZ,  gsize=KA,    start=startZ )
    if ( IO_PNETCDF ) then
       call HistoryPutAxis( 'CDX',  'Grid Cell length X', 'm', 'CX',  GRID_CDXG, gsize=IAG,   start=startZ )
       call HistoryPutAxis( 'CDY',  'Grid Cell length Y', 'm', 'CY',  GRID_CDYG, gsize=JAG,   start=startZ )
    else
       call HistoryPutAxis( 'CDX',  'Grid Cell length X', 'm', 'CX',  GRID_CDX )
       call HistoryPutAxis( 'CDY',  'Grid Cell length Y', 'm', 'CY',  GRID_CDY )
    end if
    call HistoryPutAxis( 'FDZ',  'Grid distance Z',    'm', 'FDZ', GRID_FDZ,  gsize=KA-1,  start=startZ )
    if ( IO_PNETCDF ) then
       call HistoryPutAxis( 'FDX',  'Grid distance X',    'm', 'FDX', GRID_FDXG, gsize=IAG-1, start=startZ )
       call HistoryPutAxis( 'FDY',  'Grid distance Y',    'm', 'FDY', GRID_FDYG, gsize=JAG-1, start=startZ )
    else
       call HistoryPutAxis( 'FDX',  'Grid distance X',    'm', 'FDX', GRID_FDX )
       call HistoryPutAxis( 'FDY',  'Grid distance Y',    'm', 'FDY', GRID_FDY )
    end if

    call HistoryPutAxis( 'LCZ',  'Land Grid Center Position Z', 'm', 'LCZ', GRID_LCZ, down=.true., gsize=LKMAX,   start=startZ )
    call HistoryPutAxis( 'LFZ',  'Land Grid Face Position Z',   'm', 'LFZ', GRID_LFZ, down=.true., gsize=LKMAX+1, start=startZ )
    call HistoryPutAxis( 'LCDZ', 'Land Grid Cell length Z',     'm', 'LCZ', GRID_LCDZ,             gsize=LKMAX,   start=startZ )

    call HistoryPutAxis( 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', GRID_UCZ, down=.true., gsize=UKMAX,   start=startZ )
    call HistoryPutAxis( 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', GRID_UFZ, down=.true., gsize=UKMAX+1, start=startZ )
    call HistoryPutAxis( 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', GRID_UCDZ,             gsize=UKMAX,   start=startZ )

    call HistoryPutAxis('CBFZ',  'Boundary factor Center Z', '1', 'CZ', GRID_CBFZ,  gsize=KA,  start=startZ )
    if ( IO_PNETCDF ) then
       call HistoryPutAxis('CBFX',  'Boundary factor Center X', '1', 'CX', GRID_CBFXG, gsize=IAG, start=startZ )
       call HistoryPutAxis('CBFY',  'Boundary factor Center Y', '1', 'CY', GRID_CBFYG, gsize=JAG, start=startZ )
    else
       call HistoryPutAxis('CBFX',  'Boundary factor Center X', '1', 'CX', GRID_CBFX )
       call HistoryPutAxis('CBFY',  'Boundary factor Center Y', '1', 'CY', GRID_CBFY )
    end if
    call HistoryPutAxis('FBFZ',  'Boundary factor Face Z',   '1', 'CZ', GRID_FBFZ,  gsize=KA,  start=startZ )
    if ( IO_PNETCDF ) then
       call HistoryPutAxis('FBFX',  'Boundary factor Face X',   '1', 'CX', GRID_FBFXG, gsize=IAG, start=startZ )
       call HistoryPutAxis('FBFY',  'Boundary factor Face Y',   '1', 'CY', GRID_FBFYG, gsize=JAG, start=startZ )
    else
       call HistoryPutAxis('FBFX',  'Boundary factor Face X',   '1', 'CX', GRID_FBFX )
       call HistoryPutAxis('FBFY',  'Boundary factor Face Y',   '1', 'CY', GRID_FBFY )
    end if

    ! TODO: skip 8 axes below when IO_PNETCDF is true, as all axes are now global
    call HistoryPutAxis('CXG',   'Grid Center Position X (global)', 'm', 'CXG', GRID_CXG, gsize=IAG,   start=startZ )
    call HistoryPutAxis('CYG',   'Grid Center Position Y (global)', 'm', 'CYG', GRID_CYG, gsize=JAG,   start=startZ )
    call HistoryPutAxis('FXG',   'Grid Face Position X (global)',   'm', 'FXG', GRID_FXG, gsize=IAG+1, start=startZ )
    call HistoryPutAxis('FYG',   'Grid Face Position Y (global)',   'm', 'FYG', GRID_FYG, gsize=JAG+1, start=startZ )

    call HistoryPutAxis('CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', GRID_CBFXG, gsize=IAG, start=startZ )
    call HistoryPutAxis('CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', GRID_CBFYG, gsize=JAG, start=startZ )
    call HistoryPutAxis('FBFXG', 'Boundary factor Face X (global)',   '1', 'CXG', GRID_FBFXG, gsize=IAG, start=startZ )
    call HistoryPutAxis('FBFYG', 'Boundary factor Face Y (global)',   '1', 'CYG', GRID_FBFYG, gsize=JAG, start=startZ )

    ! associate coordinates
    if ( HIST_BND ) then
       start(1) = ISGB   ! global subarray starting index
       start(2) = JSGB   ! global subarray starting index
    else
       start(1) = 1 + PRC_2Drank(PRC_myrank,1) * IMAX    ! no IHALO
       start(2) = 1 + PRC_2Drank(PRC_myrank,2) * JMAX    ! no JHALO
    end if
    start(3) = 1

    do k = 1, KMAX
       AXIS(1:im,1:jm,k) = REAL_CZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ', 'y ', 'z '/)
    call HistoryPutAssociatedCoordinates( 'height', 'height above ground level', &
                                          'm', AXIS_name(1:3), AXIS(:,:,:), start=start )

    do k = 1, KMAX
       AXIS(1:im,1:jm,k) = REAL_FZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ', 'y ', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_xyw', 'height above ground level (half level xyw)', &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:), start=start )

    do k = 1, KMAX
    do j = 1, jm
    do i = 1, im-1
       AXIS(i,j,k) = ( REAL_CZ(k+KS-1,ims+i-1,jms+j-1) + REAL_CZ(k+KS-1,ims+i,jms+j-1) ) * 0.5_RP
    enddo
    enddo
    enddo
    do k = 1, KMAX
    do j = 1, jm
       AXIS(im,j,k) = REAL_CZ(k+KS-1,ims+im-1,jms+j-1)
    enddo
    enddo
    AXIS_name(1:3) = (/'xh', 'y ', 'z '/)
    call HistoryPutAssociatedCoordinates( 'height_uyz', 'height above ground level (half level uyz)', &
                                          'm', AXIS_name(1:3), AXIS(:,:,:), start=start )

    do k = 1, KMAX
    do j = 1, jm-1
    do i = 1, im
       AXIS(i,j,k) = ( REAL_CZ(k+KS-1,ims+i-1,jms+j-1) + REAL_CZ(k+KS-1,ims+i-1,jms+j) ) * 0.5_RP
    enddo
    enddo
    enddo
    do k = 1, KMAX
    do i = 1, im
       AXIS(i,jm,k) = REAL_CZ(k+KS-1,ims+i-1,jms+jm-1)
    enddo
    enddo
    AXIS_name(1:3) = (/'x ', 'yh', 'z '/)
    call HistoryPutAssociatedCoordinates( 'height_xvz', 'height above ground level (half level xvz)', &
                                          'm', AXIS_name(1:3), AXIS(:,:,:), start=start )

    do k = 1, KMAX
    do j = 1, jm-1
    do i = 1, im-1
       AXIS(i,j,k) = ( REAL_CZ(k+KS-1,ims+i-1,jms+j-1) + REAL_CZ(k+KS-1,ims+i  ,jms+j-1) &
                     + REAL_CZ(k+KS-1,ims+i-1,jms+j  ) + REAL_CZ(k+KS-1,ims+i  ,jms+j  ) ) * 0.25_RP
    enddo
    enddo
    enddo
    do k = 1, KMAX
    do i = 1, im
       AXIS(i,jm,k) = REAL_CZ(k+KS-1,ims+i-1,jms+jm-1)
    enddo
    enddo
    do k = 1, KMAX
    do j = 1, jm
       AXIS(im,j,k) = REAL_CZ(k+KS-1,ims+im-1,jms+j-1)
    enddo
    enddo
    AXIS_name(1:3) = (/'xh', 'yh', 'z '/)
    call HistoryPutAssociatedCoordinates( 'height_uvz', 'height above ground level (half level uvz)', &
                                          'm', AXIS_name(1:3), AXIS(:,:,:), start=start )

    do k = 1, KMAX
    do j = 1, jm
    do i = 1, im-1
       AXIS(i,j,k) = ( REAL_FZ(k+KS-1,ims+i-1,jms+j-1) + REAL_FZ(k+KS-1,ims+i,jms+j-1) ) * 0.5_RP
    enddo
    enddo
    enddo
    do k = 1, KMAX
    do j = 1, jm
       AXIS(im,j,k) = REAL_FZ(k+KS-1,ims+im-1,jms+j-1)
    enddo
    enddo
    AXIS_name(1:3) = (/'xh', 'y ', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_uyw', 'height above ground level (half level uyw)', &
                                          'm', AXIS_name(1:3), AXIS(:,:,:), start=start )

    do k = 1, KMAX
    do j = 1, jm-1
    do i = 1, im
       AXIS(i,j,k) = ( REAL_FZ(k+KS-1,ims+i-1,jms+j-1) + REAL_FZ(k+KS-1,ims+i-1,jms+j) ) * 0.5_RP
    enddo
    enddo
    enddo
    do k = 1, KMAX
    do i = 1, im
       AXIS(i,jm,k) = REAL_FZ(k+KS-1,ims+i-1,jms+jm-1)
    enddo
    enddo
    AXIS_name(1:3) = (/'x ', 'yh', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_xvw', 'height above ground level (half level xvw)', &
                                          'm', AXIS_name(1:3), AXIS(:,:,:), start=start )

    do k = 1, KMAX
    do j = 1, jm-1
    do i = 1, im-1
       AXIS(i,j,k) = ( REAL_FZ(k+KS-1,ims+i-1,jms+j-1) + REAL_FZ(k+KS-1,ims+i  ,jms+j-1) &
                     + REAL_FZ(k+KS-1,ims+i-1,jms+j  ) + REAL_FZ(k+KS-1,ims+i  ,jms+j  ) ) * 0.25_RP
    enddo
    enddo
    enddo
    do k = 1, KMAX
    do i = 1, im
       AXIS(i,jm,k) = REAL_FZ(k+KS-1,ims+i-1,jms+jm-1)
    enddo
    enddo
    do k = 1, KMAX
    do j = 1, jm
       AXIS(im,j,k) = REAL_FZ(k+KS-1,ims+im-1,jms+j-1)
    enddo
    enddo
    AXIS_name(1:3) = (/'xh', 'yh', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_uvw', 'height above ground level (half level uvw)', &
                                          'm', AXIS_name(1:3), AXIS(:,:,:), start=start               )

    AXIS(1:im,1:jm,1) = REAL_LON (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lon', 'longitude',                                      &
                                          'degrees_east', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = REAL_LONX(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lon_uy', 'longitude (half level uy)',                   &
                                          'degrees_east', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = REAL_LONY(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'yh'/)
    call HistoryPutAssociatedCoordinates( 'lon_xv', 'longitude (half level xv)',                   &
                                          'degrees_east', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = REAL_LONXY(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'yh'/)
    call HistoryPutAssociatedCoordinates( 'lon_uv', 'longitude (half level uv)',                   &
                                          'degrees_east', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = REAL_LAT (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lat', 'latitude',                                        &
                                          'degrees_north', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = REAL_LATX(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lat_uy', 'latitude (half level uy)',                     &
                                          'degrees_north', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = REAL_LATY(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'yh'/)
    call HistoryPutAssociatedCoordinates( 'lat_xv', 'latitude (half level xv)',                     &
                                          'degrees_north', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = REAL_LATXY(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'yh'/)
    call HistoryPutAssociatedCoordinates( 'lat_uv', 'latitude (half level uv)',                     &
                                          'degrees_north', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = TOPO_Zsfc(ims:ime,jms:jme)
    AXIS_name(1:2) = (/'x ', 'y '/)
    call HistoryPutAssociatedCoordinates( 'topo', 'topography',                         &
                                          'm', AXIS_name(1:2), AXIS(:,:,1), start=start )

    AXIS(1:im,1:jm,1) = LANDUSE_frac_land(ims:ime,jms:jme)
    AXIS_name(1:2) = (/'x ', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lsmask', 'fraction for land-sea mask', &
                                          '0-1', AXIS_name(1:2), AXIS(:,:,1), start=start )

    return
  end subroutine HIST_put_axes

  !-----------------------------------------------------------------------------
  !> set switch
  subroutine HIST_switch( switch )
    implicit none

    logical, intent(in) :: switch
    !---------------------------------------------------------------------------

    enabled = switch

    return
  end subroutine HIST_switch

  !-----------------------------------------------------------------------------
  !> set interpolation factor for pressure coordinate
  subroutine HIST_setpres( &
       PRES,    &
       SFC_PRES )
    use scale_interpolation, only: &
       INTERP_update_pres
    implicit none

    real(RP), intent(in) :: PRES    (KA,IA,JA) ! pressure in Xi coordinate [Pa]
    real(RP), intent(in) :: SFC_PRES(   IA,JA) ! surface pressure          [Pa]
    !---------------------------------------------------------------------------

    if ( HIST_PRES_nlayer > 0 ) then
       call INTERP_update_pres( HIST_PRES_nlayer,     & ! [IN]
                                PRES         (:,:,:), & ! [IN]
                                SFC_PRES     (:,:)  , & ! [IN]
                                HIST_PRES_val(:)      ) ! [IN]
    endif

    return
  end subroutine HIST_setpres

  !-----------------------------------------------------------------------------
  !> Register/Append variable to history file
  subroutine HIST_reg( &
       itemid, &
       item,   &
       desc,   &
       unit,   &
       ndim,   &
       xdim,   &
       ydim,   &
       zdim    )
    use gtool_history, only: &
       HistoryAddVariable, &
       HistoryCheck
    use scale_time, only: &
       NOWSTEP => TIME_NOWSTEP
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    use scale_process, only: &
       PRC_myrank, &
       PRC_LOCAL_COMM_WORLD
    use MPI, only : MPI_COMM_NULL
    implicit none

    integer,          intent(out) :: itemid !< index number of the item
    character(len=*), intent(in)  :: item   !< name         of the item
    character(len=*), intent(in)  :: desc   !< description  of the item
    character(len=*), intent(in)  :: unit   !< unit         of the item
    integer,          intent(in)  :: ndim   !< dimension    of the item
    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim

    logical :: flag_half_x
    logical :: flag_half_y
    logical :: flag_half_z

    character(len=H_SHORT) :: dims(3)

    integer                :: nvariant1, nvariant2, nvariant3
    integer                :: v, id
    logical                :: atom

    integer                :: rankidx(2)
    integer                :: start(4), count(4)
    integer                :: comm

    !---------------------------------------------------------------------------

    itemid = -1

    if( .NOT. enabled ) return

    if( HIST_item_limit == 0 ) return

    do id = 1, HIST_item_count
       if ( item == HIST_item(id) ) then ! item exists
          itemid = id
          return
       endif
    enddo

    call PROF_rapstart('FILE_O_NetCDF', 2)

    ! Try to add new item

    atom = .true.

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start = 0
    count = 0

    if ( ndim == 1 ) then

       ! check half/full level for vertical
       dims(1) = "z"
       if ( present(zdim) ) then
          if ( zdim == 'half' ) then
             dims(1) = "zh"
          endif
       endif

       ! for shared-file parallel I/O, only rank 0 writes variables with only Z dimension
       start(1) = 1
       count(1) = KMAX
       if ( PRC_myrank .GT. 0 ) count(1) = 0

    elseif ( ndim == 2 ) then

       ! check half/full level for horizontal
       flag_half_x = .false.
       if ( present(xdim) ) then
          if( xdim == 'half' ) flag_half_x = .true.
       endif

       flag_half_y = .false.
       if ( present(ydim) ) then
          if( ydim == 'half' ) flag_half_y = .true.
       endif

       if    ( flag_half_x .AND. flag_half_y ) then
          dims(1) = 'lon_uv'
          dims(2) = 'lat_uv'
       elseif( flag_half_x ) then
          dims(1) = 'lon_uy'
          dims(2) = 'lat_uy'
       elseif( flag_half_y ) then
          dims(1) = 'lon_xv'
          dims(2) = 'lat_xv'
       else
          dims(1) = 'lon'
          dims(2) = 'lat'
       endif

    elseif ( ndim == 3 ) then

       ! check half/full level for vertical/horizontal
       flag_half_x = .false.

       if ( present(xdim) ) then
          if( xdim == 'half' ) flag_half_x = .true.
       endif

       flag_half_y = .false.
       if ( present(ydim) ) then
          if( ydim == 'half' ) flag_half_y = .true.
       endif

       flag_half_z = .false.
       if ( present(zdim) ) then
          if( zdim == 'half' ) flag_half_z = .true.
       endif

       if    ( flag_half_x .AND. flag_half_y ) then
          dims(1) = 'lon_uv'
          dims(2) = 'lat_uv'
          if ( flag_half_z ) then
             dims(3) = 'height_uvw'
          else
             dims(3) = 'height_uvz'
          endif
       elseif( flag_half_x ) then
          dims(1) = 'lon_uy'
          dims(2) = 'lat_uy'
          if ( flag_half_z ) then
             dims(3) = 'height_uyw'
          else
             dims(3) = 'height_uyz'
          endif
       elseif( flag_half_y ) then
          dims(1) = 'lon_xv'
          dims(2) = 'lat_xv'
          if ( flag_half_z ) then
             dims(3) = 'height_xvw'
          else
             dims(3) = 'height_xvz'
          endif
       else
          dims(1) = 'lon'
          dims(2) = 'lat'
          if ( flag_half_z ) then
             dims(3) = 'height_xyw'
          else
             dims(3) = 'height'
          endif
       endif

       ! start and count will be used by PnetCDF I/O
       start(3) = 1
       count(3) = KMAX

       if ( present(zdim) ) then
          if    ( zdim == 'land'      ) then
             dims(3) = 'lz'
             count(3) = LKMAX
             atom = .false.
          elseif( zdim == 'landhalf'  ) then
             dims(3) = 'lzh'
             count(3) = LKMAX
             atom = .false.
          elseif( zdim == 'urban'     ) then
             dims(3) = 'uz'
             count(3) = UKMAX
             atom = .false.
          elseif( zdim == 'urbanhalf' ) then
             dims(3) = 'uzh'
             count(3) = UKMAX
             atom = .false.
          endif

       endif

    endif

    if ( ndim >= 2 ) then
       ! start and count will be used by PnetCDF I/O
       if ( HIST_BND ) then
          start(1) = ISGB
          start(2) = JSGB
          count(1) = IMAXB
          count(2) = JMAXB
       else
          ! for the case the shared-file contains no halos
          start(1) = 1 + PRC_2Drank(PRC_myrank,1) * IMAX    ! no IHALO
          start(2) = 1 + PRC_2Drank(PRC_myrank,2) * JMAX    ! no JHALO
          count(1) = IMAX
          count(2) = JMAX
       end if
    end if


    if ( IO_PNETCDF ) then  ! user input parameter indicates to do PnetCDF I/O
       comm = PRC_LOCAL_COMM_WORLD
    else
       comm = MPI_COMM_NULL
    end if

    if ( atom ) then

       ! model coordinate (terrain following coordinate)
       call HistoryAddVariable( nvariant1,    & ! [OUT]
                                item,         & ! [IN]
                                dims(1:ndim), & ! [IN]
                                desc,         & ! [IN]
                                unit,         & ! [IN]
                                NOWSTEP,      & ! [IN]
                                'model',      & ! [IN]
                                start=start,  & ! [IN]
                                count=count,  & ! [IN]
                                comm=comm     ) ! [IN]

       ! absolute height coordinate
       dims(3) = 'z'
       call HistoryAddVariable( nvariant2,    & ! [OUT]
                                item,         & ! [IN]
                                dims(1:ndim), & ! [IN]
                                desc,         & ! [IN]
                                unit,         & ! [IN]
                                NOWSTEP,      & ! [IN]
                                'z',          & ! [IN]
                                start=start,  & ! [IN]
                                count=count,  & ! [IN]
                                comm=comm     ) ! [IN]

       ! pressure coordinate
       if ( HIST_PRES_nlayer > 0 ) then

          dims(3) = 'pressure'
          call HistoryAddVariable( nvariant3,    & ! [OUT]
                                   item,         & ! [IN]
                                   dims(1:ndim), & ! [IN]
                                   desc,         & ! [IN]
                                   unit,         & ! [IN]
                                   NOWSTEP,      & ! [IN]
                                   'pressure',   & ! [IN]
                                   start=start,  & ! [IN]
                                   count=count,  & ! [IN]
                                   comm=comm     ) ! [IN]

       else
          nvariant3 = 0
       endif

    else

       call HistoryAddVariable( nvariant1,    & ! [OUT]
                                item,         & ! [IN]
                                dims(1:ndim), & ! [IN]
                                desc,         & ! [IN]
                                unit,         & ! [IN]
                                NOWSTEP,      & ! [IN]
                                start=start,  & ! [IN]
                                count=count,  & ! [IN]
                                comm=comm     ) ! [IN]

       nvariant2 = 0
       nvariant3 = 0
    endif

    if ( nvariant1 + nvariant2 + nvariant3 > 0 ) then
       HIST_item_count   = HIST_item_count + 1
       itemid            = HIST_item_count
       HIST_item(itemid) = item

       do v = 1, nvariant1
          HIST_variant(itemid)                      = HIST_variant(itemid) + 1
          HIST_zcoord (itemid,HIST_variant(itemid)) = I_MODEL
       enddo

       do v = 1, nvariant2
          HIST_variant(itemid)                      = HIST_variant(itemid) + 1
          HIST_zcoord (itemid,HIST_variant(itemid)) = I_Z
       enddo

       do v = 1, nvariant3
          HIST_variant(itemid)                      = HIST_variant(itemid) + 1
          HIST_zcoord (itemid,HIST_variant(itemid)) = I_PRES
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_reg

  !-----------------------------------------------------------------------------
  !> Check time to putting data
  subroutine HIST_query( &
       itemid, &
       answer  )
    use gtool_history, only: &
       HistoryQuery
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    integer,  intent(in)  :: itemid !< name of the item
    logical,  intent(out) :: answer !< is it time to store?
    !---------------------------------------------------------------------------

    answer = .false.

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    call HistoryQuery( HIST_item(itemid), & ! [IN]
                       TIME_NOWSTEP,      & ! [IN]
                       answer             ) ! [OUT]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_query

  !-----------------------------------------------------------------------------
  !> Put 1D data to history buffer
  subroutine HIST_put_0D( &
       itemid, &
       var     )
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    integer,  intent(in)  :: itemid !< name of the item
    real(RP), intent(in)  :: var    !< value

    integer  :: n, v, id
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    id = 0
    do n = 1, itemid-1
       id = id + HIST_variant(n)
    enddo

    do v = 1, HIST_variant(itemid)
       id = id + 1

       call HistoryPut( id,           & ! [IN]
                        TIME_NOWSTEP, & ! [IN]
                        var           ) ! [IN]
    enddo

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_put_0D

  !-----------------------------------------------------------------------------
  !> Put 1D data to history buffer
  subroutine HIST_put_1D( &
       itemid, &
       var     )
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    integer,  intent(in)  :: itemid !< name of the item
    real(RP), intent(in)  :: var(:) !< value

    real(RP) :: var_trim(KMAX)

    integer  :: n, v, id
    integer  :: k
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    do k = 1, KMAX
       var_trim(k) = var(KS+k-1)
    enddo

    id = 0
    do n = 1, itemid-1
       id = id + HIST_variant(n)
    enddo

    do v = 1, HIST_variant(itemid)
       id = id + 1

       call HistoryPut( id,           & ! [IN]
                        TIME_NOWSTEP, & ! [IN]
                        var_trim(:)   ) ! [IN]
    enddo

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_put_1D

  !-----------------------------------------------------------------------------
  !> Put 2D data to history buffer
  subroutine HIST_put_2D( &
       itemid, &
       var,    &
       nohalo  )
    use gtool_file, only: &
       RMISS
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    integer,  intent(in)  :: itemid   !< name of the item
    real(RP), intent(in)  :: var(:,:) !< value
    logical,  intent(in), optional :: nohalo

    real(RP) :: var_trim(im*jm)
    logical  :: nohalo_

    integer  :: n, v, id
    integer  :: i, j
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    do j = 1, jm
    do i = 1, im
       var_trim((j-1)*im+i) = var(ims+i-1,jms+j-1)
    enddo
    enddo

    nohalo_ = .false.
    if( present(nohalo) ) nohalo_ = nohalo

    if ( nohalo_ ) then
       ! W halo
       do j = 1, jm
       do i = 1, IS-ims
          var_trim((j-1)*im+i) = RMISS
       enddo
       enddo
       ! E halo
       do j = 1, jm
       do i = IE-ims+2, ime-ims+1
          var_trim((j-1)*im+i) = RMISS
       enddo
       enddo
       ! S halo
       do j = 1, JS-jms
       do i = 1, im
          var_trim((j-1)*im+i) = RMISS
       enddo
       enddo
       ! N halo
       do j = JE-jms+2, jme-jms+1
       do i = 1, im
          var_trim((j-1)*im+i) = RMISS
       enddo
       enddo
    endif

    id = 0
    do n = 1, itemid-1
       id = id + HIST_variant(n)
    enddo

    do v = 1, HIST_variant(itemid)
       id = id + 1

       call HistoryPut( id,           & ! [IN]
                        TIME_NOWSTEP, & ! [IN]
                        var_trim(:)   ) ! [IN]
    enddo

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_put_2D

  !-----------------------------------------------------------------------------
  !> Put 3D data to history buffer
  subroutine HIST_put_3D( &
       itemid, &
       var,    &
       xdim,   &
       ydim,   &
       zdim,   &
       nohalo  )
    use gtool_file, only: &
       RMISS
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWSTEP
    use scale_interpolation, only: &
       INTERP_vertical_xi2z, &
       INTERP_vertical_xi2p, &
       INTERP_available
    implicit none

    integer,          intent(in)  :: itemid     !< name of the item
    real(RP),         intent(in)  :: var(:,:,:) !< value
    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim
    logical,          intent(in), optional :: nohalo

    character(len=H_SHORT) :: xd, yd, zd
    integer                :: isize, jsize, ksize
    integer                :: istart, jstart, kstart

    real(RP) :: var_Z(KA              ,IA,JA)
    real(RP) :: var_P(HIST_PRES_nlayer,IA,JA)

    real(RP) :: var_trim(km*im*jm)
    logical  :: nohalo_
    integer  :: s(3)

    integer  :: n, v, id
    integer  :: i, j, k, ijk

    intrinsic shape
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    xd = ''
    yd = ''
    zd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim
    if( present(zdim) ) zd = zdim

    nohalo_ = .false.
    if( present(nohalo) ) nohalo_ = nohalo

    ! select dimension
    select case ( xd )
      case ('half')
        isize  = im
        istart = ims
      case default
        isize  = im
        istart = ims
    end select

    select case ( yd )
      case ('half')
        jsize  = jm
        jstart = jms
      case default
        jsize  = jm
        jstart = jms
    end select

    select case ( zd )
      case ('land')
        ksize  = LKMAX
        kstart = LKS
      case ('landhalf')
        ksize  = LKMAX
        kstart = LKS
      case ('urban')
        ksize  = UKMAX
        kstart = UKS
      case ('urbanhalf')
        ksize  = UKMAX
        kstart = UKS
      case ('half')
        ksize  = KMAX
        kstart = KS
      case default
        ksize  = KMAX
        kstart = KS
    end select

    s(:) = shape(var)

    id = 0
    do n = 1, itemid-1
       id = id + HIST_variant(n)
    enddo

    do v = 1, HIST_variant(itemid)

       if    (       s(1)  == KA                  &
               .AND. ksize == KMAX                &
               .AND. HIST_zcoord(itemid,v) == I_Z &
               .AND. INTERP_available             ) then ! z*->z interpolation

          call PROF_rapstart('FILE_O_interp', 2)
          call INTERP_vertical_xi2z( var  (:,:,:), & ! [IN]
                                     var_Z(:,:,:)  ) ! [OUT]
          call PROF_rapend  ('FILE_O_interp', 2)

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = var_Z(kstart+k-1,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo

       elseif(       s(1)  == KA                     &
               .AND. ksize == KMAX                   &
               .AND. HIST_zcoord(itemid,v) == I_PRES ) then ! z*->p interpolation

          ksize = HIST_PRES_nlayer

          call PROF_rapstart('FILE_O_interp', 2)
          call INTERP_vertical_xi2p( HIST_PRES_nlayer, & ! [IN]
                                     var  (:,:,:),     & ! [IN]
                                     var_P(:,:,:)      ) ! [OUT]
          call PROF_rapend  ('FILE_O_interp', 2)

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = var_P(k,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo

       else ! no interpolation

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = var(kstart+k-1,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo

       endif

       if ( nohalo_ ) then
             ! W halo
          do k = 1, ksize
          do j = 1, jsize
          do i = 1, IS-istart
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
          enddo
          enddo
          enddo
          ! E halo
          do k = 1, ksize
          do j = 1, jsize
          do i = IE-istart+2, ime-istart+1
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
          enddo
          enddo
          enddo
          ! S halo
          do k = 1, ksize
          do j = 1, JS-jstart
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
          enddo
          enddo
          enddo
          ! N halo
          do k = 1, ksize
          do j = JE-jstart+2, jme-jstart+1
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
          enddo
          enddo
          enddo
       endif

       id = id + 1

       call HistoryPut( id,                           & ! [IN]
                        TIME_NOWSTEP,                 & ! [IN]
                        var_trim(1:isize*jsize*ksize) ) ! [IN]
    enddo

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_put_3D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 0D
  subroutine HIST_in_0D( &
       var,  &
       item, &
       desc, &
       unit )
    implicit none

    real(RP),         intent(in) :: var
    character(len=*), intent(in) :: item
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: unit

    integer :: itemid
    logical :: do_put
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    ! Check whether the item has been already registered
    call HIST_reg  ( itemid,  & ! [OUT]
                     item,    & ! [IN]
                     desc,    & ! [IN]
                     unit,    & ! [IN]
                     0        ) ! [IN]

    ! Check whether it is time to input the item
    call HIST_query( itemid,  & ! [IN]
                     do_put   ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid, & ! [IN]
                      var     ) ! [IN]
    endif

    return
  end subroutine HIST_in_0D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 1D
  subroutine HIST_in_1D( &
       var,  &
       item, &
       desc, &
       unit, &
       zdim  )
    implicit none

    real(RP),         intent(in)  :: var(:)
    character(len=*), intent(in)  :: item
    character(len=*), intent(in)  :: desc
    character(len=*), intent(in)  :: unit
    character(len=*), intent(in), optional :: zdim

    character(len=H_SHORT) :: zd

    integer :: itemid
    logical :: do_put
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    zd = ''
    if( present(zdim) ) zd = zdim

    ! Check whether the item has been already registered
    call HIST_reg  ( itemid,  & ! [OUT]
                     item,    & ! [IN]
                     desc,    & ! [IN]
                     unit,    & ! [IN]
                     1,       & ! [IN]
                     zdim=zd  ) ! [IN]

    ! Check whether it is time to input the item
    call HIST_query( itemid,  & ! [IN]
                     do_put   ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid, & ! [IN]
                      var     ) ! [IN]
    endif

    return
  end subroutine HIST_in_1D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 2D
  subroutine HIST_in_2D( &
       var,   &
       item,  &
       desc,  &
       unit,  &
       xdim,  &
       ydim,  &
       nohalo )
    implicit none

    real(RP),         intent(in)  :: var(:,:) !< value
    character(len=*), intent(in)  :: item     !< name        of the item
    character(len=*), intent(in)  :: desc     !< description of the item
    character(len=*), intent(in)  :: unit     !< unit        of the item
    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    logical,          intent(in), optional :: nohalo

    character(len=H_SHORT) :: xd, yd

    integer :: itemid
    logical :: do_put
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    xd = ''
    yd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim

    ! Check whether the item has been already registered
    call HIST_reg  ( itemid,  & ! [OUT]
                     item,    & ! [IN]
                     desc,    & ! [IN]
                     unit,    & ! [IN]
                     2,       & ! [IN]
                     xdim=xd, & ! [IN]
                     ydim=yd  ) ! [IN]

    ! Check whether it is time to input the item
    call HIST_query( itemid,  & ! [IN]
                     do_put   ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid,       & ! [IN]
                      var,          & ! [IN]
                      nohalo=nohalo ) ! [IN]
    endif

    return
  end subroutine HIST_in_2D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 3D
  subroutine HIST_in_3D( &
       var,   &
       item,  &
       desc,  &
       unit,  &
       xdim,  &
       ydim,  &
       zdim,  &
       nohalo )
    implicit none

    real(RP),         intent(in)  :: var(:,:,:) !< value
    character(len=*), intent(in)  :: item       !< name        of the item
    character(len=*), intent(in)  :: desc       !< description of the item
    character(len=*), intent(in)  :: unit       !< unit        of the item
    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim
    logical,          intent(in), optional :: nohalo

    character(len=H_SHORT) :: xd, yd, zd

    integer :: itemid
    logical :: do_put
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    xd = ''
    yd = ''
    zd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim
    if( present(zdim) ) zd = zdim

    ! Check whether the item has been already registered
    call HIST_reg  ( itemid,  & ! [OUT]
                     item,    & ! [IN]
                     desc,    & ! [IN]
                     unit,    & ! [IN]
                     3,       & ! [IN]
                     xdim=xd, & ! [IN]
                     ydim=yd, & ! [IN]
                     zdim=zd  ) ! [IN]

    ! Check whether it is time to input the item
    call HIST_query( itemid,  & ! [IN]
                     do_put   ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid,       & ! [IN]
                      var,          & ! [IN]
                      xdim=xd,      & ! [IN]
                      ydim=yd,      & ! [IN]
                      zdim=zd,      & ! [IN]
                      nohalo=nohalo ) ! [IN]
    endif

    return
  end subroutine HIST_in_3D

  !-----------------------------------------------------------------------------
  !> Get 1D data from file
  subroutine HIST_get_1D( &
       var,          &
       basename,     &
       varname,      &
       step,         &
       allow_missing )
    use gtool_history, only: &
       HistoryGet
    implicit none

    real(RP),         intent(out) :: var(:)     !< value
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    integer,          intent(in)  :: step       !< step number
    logical,          intent(in), optional :: allow_missing !< allow data is missing?

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:),          & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine HIST_get_1D

  !-----------------------------------------------------------------------------
  !> Get 2D data from file
  subroutine HIST_get_2D( &
       var,          &
       basename,     &
       varname,      &
       step,         &
       allow_missing )
    use gtool_history, only: &
       HistoryGet
    implicit none

    real(RP),         intent(out) :: var(:,:)   !< value
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    integer,          intent(in)  :: step       !< step number
    logical,          intent(in), optional :: allow_missing !< allow data is missing?

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:,:),        & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine HIST_get_2D

  !-----------------------------------------------------------------------------
  !> Get 3D data from file
  subroutine HIST_get_3D( &
       var,          &
       basename,     &
       varname,      &
       step,         &
       allow_missing )
    use gtool_history, only: &
       HistoryGet
    implicit none

    real(RP),         intent(out) :: var(:,:,:) !< value
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    integer,          intent(in)  :: step       !< step number
    logical,          intent(in), optional :: allow_missing !< allow data is missing?

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:,:,:),      & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine HIST_get_3D

  !-----------------------------------------------------------------------------
  !> Flush history buffer to file
  subroutine HIST_write
    use gtool_history, only: &
       HistoryWriteAll
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    call HistoryWriteAll( TIME_NOWSTEP ) ![IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_write

end module scale_history
