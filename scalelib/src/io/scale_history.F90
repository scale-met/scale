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
  public :: HIST_reg
  public :: HIST_in
  public :: HIST_put
  public :: HIST_get
  public :: HIST_write

  interface HIST_in
     module procedure HIST_in_1D
     module procedure HIST_in_2D
     module procedure HIST_in_3D
  end interface HIST_in
  interface HIST_put
     module procedure HIST_put_1D
     module procedure HIST_put_2D
     module procedure HIST_put_3D
  end interface HIST_put
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
  integer :: im, jm
  integer :: ims, ime
  integer :: jms, jme
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine HIST_setup
    use gtool_history, only: &
       HistoryInit
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank, &
       PRC_HAS_W, &
       PRC_HAS_E, &
       PRC_HAS_S, &
       PRC_HAS_N
    use scale_time, only: &
       TIME_OFFSET_YEAR
    implicit none

    character(len=H_MID) :: HISTORY_H_TITLE = 'SCALE-LES HISTORY OUTPUT' !< title of the output file
    character(len=H_SHORT) :: HISTORY_T_UNITS = 'seconds'
    character(len=H_MID)   :: HISTORY_T_SINCE = ''

    logical :: HIST_BND = .true.

    namelist / PARAM_HIST /      &
         HIST_BND

    integer :: rankidx(2)
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[HISTORY] / Categ[IO] / Origin[SCALElib]'

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


    call PROF_rapstart('FILE O NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank, 1)
    rankidx(2) = PRC_2Drank(PRC_myrank, 2)

    if ( TIME_OFFSET_YEAR > 0 ) then
       HISTORY_T_SINCE = '0000-01-01 00:00:00'
       write(HISTORY_T_SINCE(1:4),'(i4)') TIME_OFFSET_YEAR
    end if

    if ( HIST_BND ) then
       im = IMAXB
       jm = JMAXB
       ims = ISB
       ime = IEB
       jms = JSB
       jme = JEB
    else
       im = IMAX
       jm = JMAX
       ims = IS
       ime = IE
       jms = JS
       jme = JE
    end if

    call HistoryInit( HISTORY_H_TITLE,           &
                      H_SOURCE,                  &
                      H_INSTITUTE,               &
                      im*jm*KMAX,                &
                      PRC_master,                &
                      PRC_myrank,                &
                      rankidx,                   &
                      namelist_fid = IO_FID_CONF, &
                      time_units = HISTORY_T_UNITS, &
                      time_since = HISTORY_T_SINCE )

    call HIST_put_axes

    call PROF_rapend  ('FILE O NetCDF', 2)

    return
  end subroutine HIST_setup

  !-----------------------------------------------------------------------------
  !> Put axis coordinate to history file
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
    implicit none

    real(RP)         :: AXIS     (im,jm,KMAX)
    character(len=2) :: AXIS_name(3)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call HistoryPutAxis( 'x',   'X',               'm', 'x',   GRID_CX(ims:ime) )
    call HistoryPutAxis( 'y',   'Y',               'm', 'y',   GRID_CY(jms:jme) )
    call HistoryPutAxis( 'z',   'Z',               'm', 'z',   GRID_CZ(KS:KE) )
    call HistoryPutAxis( 'xh',  'X (half level)',  'm', 'xh',  GRID_FX(ims:ime) )
    call HistoryPutAxis( 'yh',  'Y (half level)',  'm', 'yh',  GRID_FY(jms:jme) )
    call HistoryPutAxis( 'zh',  'Z (half level)',  'm', 'zh',  GRID_FZ(KS:KE) )

    call HistoryPutAxis( 'lz',  'LZ',              'm', 'lz',  GRID_LCZ(LKS:LKE), down=.true. )
    call HistoryPutAxis( 'lzh', 'LZ (half level)', 'm', 'lzh', GRID_LFZ(LKS:LKE), down=.true. )

    call HistoryPutAxis( 'uz',  'UZ',              'm', 'uz',  GRID_UCZ(UKS:UKE), down=.true. )
    call HistoryPutAxis( 'uzh', 'UZ (half level)', 'm', 'uzh', GRID_UFZ(UKS:UKE), down=.true. )

    call HistoryPutAxis( 'CZ',  'Atmos Grid Center Position Z', 'm', 'CZ',  GRID_CZ )
    call HistoryPutAxis( 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  GRID_CX )
    call HistoryPutAxis( 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  GRID_CY )
    call HistoryPutAxis( 'FZ',  'Atmos Grid Face Position Z',   'm', 'FZ',  GRID_FZ )
    call HistoryPutAxis( 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  GRID_FX )
    call HistoryPutAxis( 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  GRID_FY )

    call HistoryPutAxis( 'CDZ',  'Grid Cell length Z', 'm', 'CZ',  GRID_CDZ )
    call HistoryPutAxis( 'CDX',  'Grid Cell length X', 'm', 'CX',  GRID_CDX )
    call HistoryPutAxis( 'CDY',  'Grid Cell length Y', 'm', 'CY',  GRID_CDY )
    call HistoryPutAxis( 'FDZ',  'Grid distance Z',    'm', 'FDZ', GRID_FDZ )
    call HistoryPutAxis( 'FDX',  'Grid distance X',    'm', 'FDX', GRID_FDX )
    call HistoryPutAxis( 'FDY',  'Grid distance Y',    'm', 'FDY', GRID_FDY )

    call HistoryPutAxis( 'LCZ',  'Land Grid Center Position Z', 'm', 'LCZ', GRID_LCZ, down=.true. )
    call HistoryPutAxis( 'LFZ',  'Land Grid Face Position Z',   'm', 'LFZ', GRID_LFZ, down=.true. )
    call HistoryPutAxis( 'LCDZ', 'Land Grid Cell length Z',     'm', 'LCZ', GRID_LCDZ )

    call HistoryPutAxis( 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', GRID_UCZ, down=.true.  )
    call HistoryPutAxis( 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', GRID_UFZ , down=.true. )
    call HistoryPutAxis( 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', GRID_UCDZ )

    call HistoryPutAxis('CBFZ',  'Boundary factor Center Z', '1', 'CZ', GRID_CBFZ)
    call HistoryPutAxis('CBFX',  'Boundary factor Center X', '1', 'CX', GRID_CBFX)
    call HistoryPutAxis('CBFY',  'Boundary factor Center Y', '1', 'CY', GRID_CBFY)
    call HistoryPutAxis('FBFZ',  'Boundary factor Face Z',   '1', 'CZ', GRID_FBFZ)
    call HistoryPutAxis('FBFX',  'Boundary factor Face X',   '1', 'CX', GRID_FBFX)
    call HistoryPutAxis('FBFY',  'Boundary factor Face Y',   '1', 'CY', GRID_FBFY)

    call HistoryPutAxis('CXG',   'Grid Center Position X (global)', 'm', 'CXG', GRID_CXG)
    call HistoryPutAxis('CYG',   'Grid Center Position Y (global)', 'm', 'CYG', GRID_CYG)
    call HistoryPutAxis('FXG',   'Grid Face Position X (global)',   'm', 'FXG', GRID_FXG)
    call HistoryPutAxis('FYG',   'Grid Face Position Y (global)',   'm', 'FYG', GRID_FYG)

    call HistoryPutAxis('CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', GRID_CBFXG)
    call HistoryPutAxis('CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', GRID_CBFYG)
    call HistoryPutAxis('FBFXG', 'Boundary factor Face X (global)',   '1', 'CXG', GRID_FBFXG)
    call HistoryPutAxis('FBFYG', 'Boundary factor Face Y (global)',   '1', 'CYG', GRID_FBFYG)

    ! associate coordinates
    do k = 1, KMAX
       AXIS(1:im,1:jm,k) = REAL_CZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ','y ', 'z '/)
    call HistoryPutAssociatedCoordinates( 'height' , 'height'             ,     &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:) )

    do k = 1, KMAX
       AXIS(1:im,1:jm,k) = REAL_FZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ','y ', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_xyw' , 'height (half level xyw)'             ,     &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:) )

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
    AXIS_name(1:3) = (/'xh','y ', 'z '/)
    call HistoryPutAssociatedCoordinates( 'height_uyz' , 'height (half level uyz)'             ,     &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:) )

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
    AXIS_name(1:3) = (/'x ','yh', 'z '/)
    call HistoryPutAssociatedCoordinates( 'height_xvz' , 'height (half level xvz)'             ,     &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:) )

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
    AXIS_name(1:3) = (/'xh','yh', 'z '/)
    call HistoryPutAssociatedCoordinates( 'height_uvz' , 'height (half level uvz)'             ,     &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:) )

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
    AXIS_name(1:3) = (/'xh','y ', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_uyw' , 'height (half level uyw)'             ,     &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:) )

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
    AXIS_name(1:3) = (/'x ','yh', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_xvw' , 'height (half level xvw)'             ,     &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:) )

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
    AXIS_name(1:3) = (/'xh','yh', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_uvw' , 'height (half level uvw)'             ,     &
                                          'm' , AXIS_name(1:3), AXIS(:,:,:) )

    AXIS(1:im,1:jm,1) = REAL_LON (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ','y '/)
    call HistoryPutAssociatedCoordinates( 'lon' , 'longitude'             ,     &
                                          'degrees_east' , AXIS_name(1:2), AXIS(:,:,1) )

    AXIS(1:im,1:jm,1) = REAL_LONX(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh','y '/)
    call HistoryPutAssociatedCoordinates( 'lon_uy', 'longitude (half level uy)',     &
                                          'degrees_east' , AXIS_name(1:2), AXIS(:,:,1) )
    AXIS(1:im,1:jm,1) = REAL_LONY(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ','yh'/)
    call HistoryPutAssociatedCoordinates( 'lon_xv', 'longitude (half level xv)',     &
                                          'degrees_east' , AXIS_name(1:2), AXIS(:,:,1) )
    AXIS(1:im,1:jm,1) = REAL_LONXY(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh','yh'/)
    call HistoryPutAssociatedCoordinates( 'lon_uv', 'longitude (half level uv)',     &
                                          'degrees_east' , AXIS_name(1:2), AXIS(:,:,1) )

    AXIS(1:im,1:jm,1) = REAL_LAT (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ','y '/)
    call HistoryPutAssociatedCoordinates( 'lat' , 'latitude'              ,     &
                                          'degrees_north', AXIS_name(1:2), AXIS(:,:,1) )

    AXIS(1:im,1:jm,1) = REAL_LATX(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh','y '/)
    call HistoryPutAssociatedCoordinates( 'lat_uy', 'latitude (half level uy)' ,     &
                                          'degrees_north', AXIS_name(1:2), AXIS(:,:,1) )
    AXIS(1:im,1:jm,1) = REAL_LATY(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ','yh'/)
    call HistoryPutAssociatedCoordinates( 'lat_xv', 'latitude (half level xv)' ,     &
                                          'degrees_north', AXIS_name(1:2), AXIS(:,:,1) )
    AXIS(1:im,1:jm,1) = REAL_LATXY(ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh','yh'/)
    call HistoryPutAssociatedCoordinates( 'lat_uv', 'latitude (half level uv)' ,     &
                                          'degrees_north', AXIS_name(1:2), AXIS(:,:,1) )

    return
  end subroutine HIST_put_axes

  !-----------------------------------------------------------------------------
  !> Register/Append variable to history file
  subroutine HIST_reg( &
       itemid,  &
       zinterp, &
       item,    &
       desc,    &
       unit,    &
       ndim,    &
       xdim,    &
       ydim,    &
       zdim     )
    use gtool_history, only: &
       HistoryAddVariable
    use scale_time, only: &
       TIME_STARTSEC
    implicit none

    integer,          intent(out) :: itemid  !< index number of the item
    logical,          intent(out) :: zinterp !< z* -> z flag of the item
    character(len=*), intent(in)  :: item    !< name         of the item
    character(len=*), intent(in)  :: desc    !< description  of the item
    character(len=*), intent(in)  :: unit    !< unit         of the item
    integer,          intent(in)  :: ndim    !< dimension    of the item

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim

    logical :: xh
    logical :: yh

    logical :: existed

    character(len=16) :: dims(3)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF', 2)

    if ( ndim == 1 ) then

       dims(1) = "z"
       if ( present(zdim) ) then
          if ( zdim=='half' ) then
             dims(1) = "zh"
          end if
       end if

    else

       xh = .false.
       if ( present(xdim) ) then
          if ( xdim=='half' ) xh = .true.
       end if
       yh = .false.
       if ( present(ydim) ) then
          if ( ydim=='half' ) yh = .true.
       end if

       if ( xh .and. yh ) then
          dims(1) = 'lon_uv'
          dims(2) = 'lat_uv'
          dims(3) = 'height_uvz'
       else if ( xh ) then
          dims(1) = 'lon_uy'
          dims(2) = 'lat_uy'
          dims(3) = 'height_uyz'
       else if ( yh ) then
          dims(1) = 'lon_xv'
          dims(2) = 'lat_xv'
          dims(3) = 'height_xvz'
       else
          dims(1) = 'lon'
          dims(2) = 'lat'
          dims(3) = 'height'
       end if

       if ( present(zdim) ) then
          if ( zdim=='land' ) then
             dims(3) = 'lz'
          else if ( zdim=='landhalf' ) then
             dims(3) = 'lzh'
          else if ( zdim=='urban' ) then
             dims(3) = 'uz'
          else if ( zdim=='urbanhalf' ) then
             dims(3) = 'uzh'
          else if ( zdim=='half' ) then
             if ( xh .and. yh ) then
                dims(3) = 'height_uvw'
             else if ( xh ) then
                dims(3) = 'height_uyw'
             else if ( yh ) then
                dims(3) = 'height_xvw'
             else
                dims(3) = 'height_xyw'
             end if
          end if
       endif

    end if

    call HistoryAddVariable( item,          & ! [IN]
                             dims(1:ndim),  & ! [IN]
                             desc,          & ! [IN]
                             unit,          & ! [IN]
                             TIME_STARTSEC, & ! [IN]
                             itemid,        & ! [OUT]
                             zinterp,       & ! [OUT]
                             existed        ) ! [OUT]

    call PROF_rapend  ('FILE O NetCDF', 2)

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
       TIME_NOWDAYSEC, &
       TIME_DTSEC
    implicit none

    integer, intent(in)  :: itemid !< index number of the item
    logical, intent(out) :: answer !< is it time to store?

    !---------------------------------------------------------------------------

    answer = .false.

    if ( itemid < 0 ) return

    call PROF_rapstart('FILE O NetCDF', 2)

    call HistoryQuery(itemid, TIME_NOWDAYSEC+TIME_DTSEC, answer)

    call PROF_rapend  ('FILE O NetCDF', 2)

    return
  end subroutine HIST_query

  !-----------------------------------------------------------------------------
  !> Put 1D data to history buffer
  subroutine HIST_put_1D( &
       itemid, &
       var     )
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWDAYSEC, &
       TIME_DTSEC
    implicit none

    integer,  intent(in) :: itemid !< index number of the item
    real(RP), intent(in) :: var(:) !< value

    real(RP) :: var2(KMAX)
    integer  :: k
    !---------------------------------------------------------------------------

    if ( itemid < 0 ) return

    call PROF_rapstart('FILE O NetCDF', 2)

    do k = 1, KMAX
       var2(k) = var(KS+k-1)
    enddo

    call HistoryPut(itemid, TIME_NOWDAYSEC+TIME_DTSEC, var2)

    call PROF_rapend  ('FILE O NetCDF', 2)

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
       TIME_NOWDAYSEC, &
       TIME_DTSEC
    implicit none

    integer,  intent(in) :: itemid   !< index number of the item
    real(RP), intent(in) :: var(:,:) !< value

    logical,  intent(in), optional :: nohalo

    real(RP) :: var2(im*jm)
    integer  :: i, j
    logical :: nohalo_
    !---------------------------------------------------------------------------

    if ( itemid < 0 ) return

    call PROF_rapstart('FILE O NetCDF', 2)

    do j = 1, jm
    do i = 1, im
       var2(i + (j-1)*im) = var(ims+i-1,jms+j-1)
    enddo
    enddo

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    if ( nohalo_ ) then
       ! W halo
       do j = 1, jm
       do i = 1, IS-ims
          var2(i + (j-1)*im) = RMISS
       enddo
       enddo
       ! E halo
       do j = 1, jm
       do i = IE-ims+2, ime-ims+1
          var2(i + (j-1)*im) = RMISS
       enddo
       enddo
       ! S halo
       do j = 1, JS-jms
       do i = 1, im
          var2(i + (j-1)*im) = RMISS
       enddo
       enddo
       ! N halo
       do j = JE-jms+2, jme-jms+1
       do i = 1, im
          var2(i + (j-1)*im) = RMISS
       enddo
       enddo
    end if

    call HistoryPut(itemid, TIME_NOWDAYSEC+TIME_DTSEC, var2)

    call PROF_rapend  ('FILE O NetCDF', 2)

    return
  end subroutine HIST_put_2D

  !-----------------------------------------------------------------------------
  !> Put 3D data to history buffer
  subroutine HIST_put_3D( &
       itemid,  &
       var,     &
       zinterp, &
       xdim,    &
       ydim,    &
       zdim,    &
       nohalo   )
    use gtool_file, only: &
       RMISS
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWDAYSEC, &
       TIME_DTSEC
    use scale_interpolation, only: &
       INTERP_vertical_xi2z, &
       INTERP_available
    implicit none

    integer,  intent(in) :: itemid     !< index number of the item
    real(RP), intent(in) :: var(:,:,:) !< value
    logical,  intent(in) :: zinterp    !< vertical interpolation?

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim
    logical,          intent(in), optional :: nohalo

    intrinsic shape
    integer :: s(3)

    character(len=16) :: xd, yd, zd

    real(RP), allocatable :: var_Z(:,:,:)
    real(RP), allocatable :: var2 (:)

    integer :: i, j, k
    integer :: isize, jsize, ksize
    integer :: iall, jall, kall
    integer :: istart, jstart, kstart

    logical :: nohalo_
    !---------------------------------------------------------------------------

    if ( itemid < 0 ) return

    call PROF_rapstart('FILE O NetCDF', 2)

    xd = ''
    yd = ''
    zd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim
    if( present(zdim) ) zd = zdim

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo


    ! select dimension
    select case ( xd )
      case default
        isize  = im
        iall   = IA
        istart = ims
    end select
    select case ( yd )
      case default
        jsize  = jm
        jall   = JA
        jstart = jms
    end select
    select case ( zd )
      case ('land')
        ksize  = LKMAX
        kall   = LKE
        kstart = LKS
      case ('urban')
        ksize  = UKMAX
        kall   = UKE
        kstart = UKS
      case default
        ksize  = KMAX
        kall   = KA
        Kstart = KS
    end select

    allocate( var_Z(kall,iall,jall)    )
    allocate( var2 (ksize*isize*jsize) )

    s = shape(var)
    if ( s(1) == 1 ) then

       do j = 1, jsize
       do i = 1, isize
          var2(i + (j-1)*isize) = var(1,istart+i-1,jstart+j-1)
       enddo
       enddo

       if ( nohalo_ ) then
          ! W halo
          do j = 1, jsize
          do i = 1, IS-ims
             var2(i + (j-1)*isize) = RMISS
          enddo
          enddo
          ! E halo
          do j = 1, jsize
          do i = IE-ims+2, ime-ims+1
             var2(i + (j-1)*isize) = RMISS
          enddo
          enddo
          ! S halo
          do j = 1, JS-jms
          do i = 1, isize
             var2(i + (j-1)*isize) = RMISS
          enddo
          enddo
          ! N halo
          do j = JE-jms+2, jme-jms+1
          do i = 1, isize
             var2(i + (j-1)*isize) = RMISS
          enddo
          enddo
       end if

       call HistoryPut(itemid, TIME_NOWDAYSEC+TIME_DTSEC, var2(1:isize*jsize))

    else
       if (       ksize == KMAX    &
            .AND. zinterp          &
            .AND. INTERP_available ) then
          call PROF_rapstart('FILE O Interpolation', 2)
          call INTERP_vertical_xi2z( var  (:,:,:), & ! [IN]
                                     var_Z(:,:,:)  ) ! [OUT]
          call PROF_rapend  ('FILE O Interpolation', 2)

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var2(i + (j-1)*isize + (k-1)*jsize*isize) = var_Z(kstart+k-1,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo
       else
          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var2(i + (j-1)*isize + (k-1)*jsize*isize) = var(kstart+k-1,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo
       endif

       if ( nohalo_ ) then
          ! W halo
          do k = 1, ksize
          do j = 1, jsize
          do i = 1, IS-ims
             var2(i + (j-1)*isize + (k-1)*jsize*isize) = RMISS
          enddo
          enddo
          enddo
          ! E halo
          do k = 1, ksize
          do j = 1, jsize
          do i = IE-ims+2, ime-ims+1
             var2(i + (j-1)*isize + (k-1)*jsize*isize) = RMISS
          enddo
          enddo
          enddo
          ! S halo
          do k = 1, ksize
          do j = 1, JS-jms
          do i = 1, isize
             var2(i + (j-1)*isize + (k-1)*jsize*isize) = RMISS
          enddo
          enddo
          enddo
          ! N halo
          do k = 1, ksize
          do j = JE-jms+2, jme-jms+1
          do i = 1, isize
             var2(i + (j-1)*isize + (k-1)*jsize*isize) = RMISS
          enddo
          enddo
          enddo
       end if

       call HistoryPut(itemid, TIME_NOWDAYSEC+TIME_DTSEC, var2)

    endif

    call PROF_rapend  ('FILE O NetCDF', 2)

    return
  end subroutine HIST_put_3D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 1D
  subroutine HIST_in_1D( &
       var,  &
       item, &
       desc, &
       unit, &
       zdim  )
    implicit none

    real(RP),         intent(in) :: var(:)
    character(len=*), intent(in) :: item
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: unit

    character(len=*), intent(in), optional :: zdim

    character(len=16) :: zd
    integer           :: itemid
    logical           :: zinterp
    logical           :: do_put
    !---------------------------------------------------------------------------

    zd = ''
    if( present(zdim) ) zd = zdim

    call HIST_reg( itemid,              & ! [OUT]
                   zinterp,             & ! [OUT]
                   item, desc, unit, 1, & ! [IN]
                   zdim = zd            ) ! [IN]

    call HIST_query( itemid, & ! [IN]
                     do_put  ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid, var ) ! [IN]
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

    real(RP),         intent(in) :: var(:,:) !< value
    character(len=*), intent(in) :: item     !< name        of the item
    character(len=*), intent(in) :: desc     !< description of the item
    character(len=*), intent(in) :: unit     !< unit        of the item

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    logical,          intent(in), optional :: nohalo

    character(len=16) :: xd, yd
    integer           :: itemid
    logical           :: zinterp
    logical           :: do_put
    !---------------------------------------------------------------------------

    xd = ''
    yd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim

    call HIST_reg( itemid,              & ! [OUT]
                   zinterp,             & ! [OUT]
                   item, desc, unit, 2, & ! [IN]
                   xdim = xd, ydim = yd ) ! [IN]

    call HIST_query( itemid, & ! [IN]
                     do_put  ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid, var, nohalo = nohalo ) ! [IN]
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

    real(RP),         intent(in) :: var(:,:,:) !< value
    character(len=*), intent(in) :: item       !< name        of the item
    character(len=*), intent(in) :: desc       !< description of the item
    character(len=*), intent(in) :: unit       !< unit        of the item

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim
    logical,          intent(in), optional :: nohalo

    character(len=16) :: xd, yd, zd
    integer           :: itemid
    logical           :: zinterp
    logical           :: do_put
    !---------------------------------------------------------------------------

    xd = ''
    yd = ''
    zd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim
    if( present(zdim) ) zd = zdim

    call HIST_reg( itemid,                       & ! [OUT]
                   zinterp,                      & ! [OUT]
                   item, desc, unit, 3,          & ! [IN]
                   xdim = xd, ydim = yd, zdim=zd ) ! [IN]

    call HIST_query( itemid, & ! [IN]
                     do_put  ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid, var, zinterp,          & ! [IN]
                      xdim = xd, ydim = yd, zdim=zd, & ! [IN]
                      nohalo = nohalo                ) ! [IN]
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

    real(RP),         intent(out) :: var(:)   !< value
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    integer,          intent(in)  :: step     !< step number

    logical,          intent(in), optional :: allow_missing

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE I NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:),          & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE I NetCDF', 2)

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

    real(RP),         intent(out) :: var(:,:) !< value
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    integer,          intent(in)  :: step     !< step number

    logical,          intent(in), optional :: allow_missing

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE I NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:,:),        & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE I NetCDF', 2)

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

    logical,          intent(in), optional :: allow_missing

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE I NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:,:,:),      & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE I NetCDF', 2)

    return
  end subroutine HIST_get_3D

  !-----------------------------------------------------------------------------
  !> Flush history buffer to file
  subroutine HIST_write
    use gtool_history, only: &
       HistoryWriteAll
    use scale_time, only: &
       TIME_NOWDAYSEC
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF', 2)

    call HistoryWriteAll( TIME_NOWDAYSEC ) ![IN]

    call PROF_rapend  ('FILE O NetCDF', 2)

    return
  end subroutine HIST_write

end module scale_history
!-------------------------------------------------------------------------------
