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
module mod_history
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use dc_types, only: &
     DP
  use gtool_history, only: &
     HistoryInit, &
     HistoryAddVariable, &
     HistoryPutAxis, &
     HistoryPutAdditionalAxis, &
     HistoryPut, &
     HistoryGet
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
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
# include "scale-les.h"
  include "inc_precision.h"
  include "inc_index.h"

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
  public :: HIST_put_axes

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_SYSCHR), private :: HISTORY_H_TITLE     = 'SCALE3 HISTORY OUTPUT' !< for header
  character(len=IO_SYSCHR), private :: HISTORY_H_SOURCE    = 'SCALE-LES ver. VERSION'//VERSION       !< for header
  character(len=IO_SYSCHR), private :: HISTORY_H_INSTITUTE = 'AICS/RIKEN'            !< for header

  character(len=1), parameter :: HISTORY_dim_name (3) = (/'x','y','z'/)              !< for axis property
  integer,          parameter :: HISTORY_dim_size (3) = (/IMAX,JMAX,KMAX/)           !< for axis property
  character(len=1), parameter :: HISTORY_dim_desc (3) = (/'X','Y','Z'/)              !< for axis property
  character(len=1), parameter :: HISTORY_dim_unit (3) = (/'m','m','m'/)              !< for axis property
  character(len=5), parameter :: HISTORY_dim_dtype(3) = (/'REAL4','REAL4','REAL4'/)  !< for axis property

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine HIST_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank
    implicit none

    integer :: rankidx(2)
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    rankidx(1) = PRC_2Drank(PRC_myrank, 1)
    rankidx(2) = PRC_2Drank(PRC_myrank, 2)

    call HistoryInit( HISTORY_H_TITLE,           &
                      HISTORY_H_SOURCE,          &
                      HISTORY_H_INSTITUTE,       &
                      PRC_master,                &
                      PRC_myrank,                &
                      rankidx,                   &
                      HISTORY_dim_name,          &
                      HISTORY_dim_size,          &
                      HISTORY_dim_desc,          &
                      HISTORY_dim_unit,          &
                      HISTORY_dim_dtype,         &
                      namelist_fid = IO_FID_CONF )

    call TIME_rapend  ('FILE O NetCDF')

    return
  end subroutine HIST_setup

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
    implicit none

    integer,          intent(out) :: itemid !< index number of the item
    character(len=*), intent(in)  :: item   !< name         of the item
    character(len=*), intent(in)  :: desc   !< description  of the item
    character(len=*), intent(in)  :: unit   !< unit         of the item
    integer,          intent(in)  :: ndim   !< dimension    of the item

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim

    logical :: existed

    character(len=2) :: dims(3)
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    dims(1) = 'x'
    dims(2) = 'y'
    dims(3) = 'z'
    if ( present(xdim) ) then
       if ( xdim=='half' ) dims(1) = 'xh'
    endif
    if ( present(ydim) ) then
       if ( ydim=='half' ) dims(2) = 'yh'
    endif
    if ( present(zdim) ) then
       if ( zdim=='half' ) dims(3) = 'zh'
    endif

    call HistoryAddVariable( item,             & ! [IN]
                             dims(1:ndim),     & ! [IN]
                             desc,             & ! [IN]
                             unit,             & ! [IN]
                             itemid  = itemid, & ! [OUT]
                             existed = existed ) ! [OUT]

    if ( .NOT. existed ) then
       call HIST_put_axes
    endif

    call TIME_rapend  ('FILE O NetCDF')

    return
  end subroutine HIST_reg

  !-----------------------------------------------------------------------------
  !> Put 1D data to history buffer
  subroutine HIST_put_1D( &
      itemid, &
      var,    &
      dt      )
    implicit none

    integer,  intent(in) :: itemid !< index number of the item
    real(RP), intent(in) :: var(:) !< value
    real(DP), intent(in) :: dt     !< delta t [sec]

    real(RP) :: var2(KMAX)
    integer  :: k
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    do k = 1, KMAX
       var2(k) = var(KS+k-1)
    enddo
    call HistoryPut(itemid, var2, dt)

    call TIME_rapend  ('FILE O NetCDF')

    return
  end subroutine HIST_put_1D

  !-----------------------------------------------------------------------------
  !> Put 2D data to history buffer
  subroutine HIST_put_2D( &
      itemid, &
      var,    &
      dt      )
    implicit none

    integer,  intent(in) :: itemid   !< index number of the item
    real(RP), intent(in) :: var(:,:) !< value
    real(DP), intent(in) :: dt       !< delta t [sec]

    real(RP) :: var2(IMAX*JMAX)
    integer  :: i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    do j = 1, JMAX
    do i = 1, IMAX
       var2(i + (j-1)*IMAX) = var(IS+i-1,JS+j-1)
    enddo
    enddo
    call HistoryPut(itemid, var2, dt)

    call TIME_rapend  ('FILE O NetCDF')

    return
  end subroutine HIST_put_2D

  !-----------------------------------------------------------------------------
  !> Put 3D data to history buffer
  subroutine HIST_put_3D( &
      itemid, &
      var,    &
      dt      )
    implicit none

    integer,  intent(in) :: itemid     !< index number of the item
    real(RP), intent(in) :: var(:,:,:) !< value
    real(DP), intent(in) :: dt         !< delta t [sec]

    intrinsic shape
    integer :: s(3)

    real(RP) :: var2(IMAX*JMAX*KMAX)
    integer  :: i, j, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    s = shape(var)
    if ( s(1) == 1 ) then

       do j = 1, JMAX
       do i = 1, IMAX
          var2(i + (j-1)*IMAX) = var(1,IS+i-1,JS+j-1)
       enddo
       enddo
       call HistoryPut(itemid, var2(1:IMAX*JMAX), dt)

    else

       do k = 1, KMAX
       do j = 1, JMAX
       do i = 1, IMAX
          var2(i + (j-1)*IMAX + (k-1)*JMAX*IMAX) = var(KS+k-1,IS+i-1,JS+j-1)
       enddo
       enddo
       enddo
       call HistoryPut(itemid, var2, dt)

    endif

    call TIME_rapend  ('FILE O NetCDF')

    return
  end subroutine HIST_put_3D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 1D
  subroutine HIST_in_1D( &
       var,  &
       item, &
       desc, &
       unit, &
       dt,   &
       zdim  )
    implicit none

    real(RP),         intent(in) :: var(:) 
    character(len=*), intent(in) :: item   
    character(len=*), intent(in) :: desc   
    character(len=*), intent(in) :: unit   
    real(DP),         intent(in) :: dt     

    character(len=*), intent(in), optional :: zdim

    character(len=4) :: zd
    integer          :: itemid
    !---------------------------------------------------------------------------

    zd = ''
    if( present(zdim) ) zd = zdim

    call HIST_reg( itemid,              & ! [OUT]
                   item, desc, unit, 1, & ! [IN]
                   zdim = zd            ) ! [IN]

    call HIST_put( itemid, var, dt ) ! [IN]

    return
  end subroutine HIST_in_1D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 2D
  subroutine HIST_in_2D( &
       var,  &
       item, &
       desc, &
       unit, &
       dt,   &
       xdim, &
       ydim  )
    implicit none

    real(RP),         intent(in) :: var(:,:) !< value
    character(len=*), intent(in) :: item     !< name        of the item
    character(len=*), intent(in) :: desc     !< description of the item
    character(len=*), intent(in) :: unit     !< unit        of the item
    real(DP),         intent(in) :: dt       !< delta t [sec]

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim

    character(len=4) :: xd, yd
    integer          :: itemid
    !---------------------------------------------------------------------------

    xd = ''
    yd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim

    call HIST_reg( itemid,              & ! [OUT]
                   item, desc, unit, 2, & ! [IN]
                   xdim = xd, ydim = yd ) ! [IN]

    call HIST_put( itemid, var, dt ) ! [IN]

    return
  end subroutine HIST_in_2D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 3D
  subroutine HIST_in_3D( &
       var,  &
       item, &
       desc, &
       unit, &
       dt,   &
       xdim, &
       ydim, &
       zdim  )
    implicit none

    real(RP),         intent(in) :: var(:,:,:) !< value
    character(len=*), intent(in) :: item       !< name        of the item
    character(len=*), intent(in) :: desc       !< description of the item
    character(len=*), intent(in) :: unit       !< unit        of the item
    real(DP),         intent(in) :: dt         !< delta t [sec]

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim

    character(len=4) :: xd, yd, zd
    integer          :: itemid
    !---------------------------------------------------------------------------

    xd = ''
    yd = ''
    zd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim
    if( present(zdim) ) zd = zdim

    call HIST_reg( itemid,                       & ! [OUT]
                   item, desc, unit, 3,          & ! [IN]
                   xdim = xd, ydim = yd, zdim=zd ) ! [IN]

    call HIST_put( itemid, var, dt ) ! [IN]

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
    implicit none

    real(RP),         intent(out) :: var(:)   !< value
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    integer,          intent(in)  :: step     !< step number

    logical,          intent(in), optional :: allow_missing

    logical :: am
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE I NetCDF')

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:),          & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call TIME_rapend  ('FILE I NetCDF')

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
    implicit none

    real(RP),         intent(out) :: var(:,:) !< value
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    integer,          intent(in)  :: step     !< step number

    logical,          intent(in), optional :: allow_missing

    logical :: am
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE I NetCDF')

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:,:),        & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call TIME_rapend  ('FILE I NetCDF')

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
    implicit none

    real(RP),         intent(out) :: var(:,:,:) !< value
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    integer,          intent(in)  :: step       !< step number

    logical,          intent(in), optional :: allow_missing

    logical :: am
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE I NetCDF')

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:,:,:),      & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call TIME_rapend  ('FILE I NetCDF')

    return
  end subroutine HIST_get_3D

  !-----------------------------------------------------------------------------
  !> Flush history buffer to file
  subroutine HIST_write
    use gtool_history, only: &
         HistoryWriteAll
    use mod_time, only: &
         TIME_NOWDAYSEC
    implicit none
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    call HistoryWriteAll( TIME_NOWDAYSEC ) ![IN]

    call TIME_rapend  ('FILE O NetCDF')

    return
  end subroutine HIST_write

  !-----------------------------------------------------------------------------
  !> Put axis coordinate to history file
  subroutine HIST_put_axes
    use mod_grid, only: &
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
    use gtool_history, only: &
       HistoryPutAxis
    implicit none
    !---------------------------------------------------------------------------

    call HistoryPutAxis('x', GRID_CX(IS:IE))
    call HistoryPutAxis('y', GRID_CY(JS:JE))
    call HistoryPutAxis('z', GRID_CZ(KS:KE))

    call HistoryPutAdditionalAxis('xh',    'X (half level)', 'm', 'xh', GRID_FX(IS:IE))
    call HistoryPutAdditionalAxis('yh',    'Y (half level)', 'm', 'yh', GRID_FY(JS:JE))
    call HistoryPutAdditionalAxis('zh',    'Z (half level)', 'm', 'zh', GRID_FZ(KS:KE))

    call HistoryPutAdditionalAxis('CZ',    'Grid Center Position Z', 'm', 'CZ', GRID_CZ)
    call HistoryPutAdditionalAxis('CX',    'Grid Center Position X', 'm', 'CX', GRID_CX)
    call HistoryPutAdditionalAxis('CY',    'Grid Center Position Y', 'm', 'CY', GRID_CY)
    call HistoryPutAdditionalAxis('FZ',    'Grid Face Position Z',   'm', 'FZ', GRID_FZ)
    call HistoryPutAdditionalAxis('FX',    'Grid Face Position X',   'm', 'FX', GRID_FX)
    call HistoryPutAdditionalAxis('FY',    'Grid Face Position Y',   'm', 'FY', GRID_FY)

    call HistoryPutAdditionalAxis('CDZ',   'Grid Cell length Z', 'm', 'CZ', GRID_CDZ)
    call HistoryPutAdditionalAxis('CDX',   'Grid Cell length X', 'm', 'CX', GRID_CDX)
    call HistoryPutAdditionalAxis('CDY',   'Grid Cell length Y', 'm', 'CY', GRID_CDY)
    call HistoryPutAdditionalAxis('FDZ',   'Grid distance Z',    'm', 'FDZ', GRID_FDZ)
    call HistoryPutAdditionalAxis('FDX',   'Grid distance X',    'm', 'FDX', GRID_FDX)
    call HistoryPutAdditionalAxis('FDY',   'Grid distance Y',    'm', 'FDY', GRID_FDY)

    call HistoryPutAdditionalAxis('CBFZ',  'Boundary factor Center Z', '1', 'CZ', GRID_CBFZ)
    call HistoryPutAdditionalAxis('CBFX',  'Boundary factor Center X', '1', 'CX', GRID_CBFX)
    call HistoryPutAdditionalAxis('CBFY',  'Boundary factor Center Y', '1', 'CY', GRID_CBFY)
    call HistoryPutAdditionalAxis('FBFZ',  'Boundary factor Face Z',   '1', 'CZ', GRID_FBFZ)
    call HistoryPutAdditionalAxis('FBFX',  'Boundary factor Face X',   '1', 'CX', GRID_FBFX)
    call HistoryPutAdditionalAxis('FBFY',  'Boundary factor Face Y',   '1', 'CY', GRID_FBFY)

    call HistoryPutAdditionalAxis('CXG',   'Grid Center Position X (global)', 'm', 'CXG', GRID_CXG)
    call HistoryPutAdditionalAxis('CYG',   'Grid Center Position Y (global)', 'm', 'CYG', GRID_CYG)
    call HistoryPutAdditionalAxis('FXG',   'Grid Face Position X (global)',   'm', 'FXG', GRID_FXG)
    call HistoryPutAdditionalAxis('FYG',   'Grid Face Position Y (global)',   'm', 'FYG', GRID_FYG)

    call HistoryPutAdditionalAxis('CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', GRID_CBFXG)
    call HistoryPutAdditionalAxis('CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', GRID_CBFYG)
    call HistoryPutAdditionalAxis('FBFXG', 'Boundary factor Face X (global)',   '1', 'CXG', GRID_FBFXG)
    call HistoryPutAdditionalAxis('FBFYG', 'Boundary factor Face Y (global)',   '1', 'CYG', GRID_FBFYG)

    return
  end subroutine HIST_put_axes

end module mod_history
!-------------------------------------------------------------------------------
