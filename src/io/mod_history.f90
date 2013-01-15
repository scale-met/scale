!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!> module History
!!
!! @par Description
!!          History output module
!!
!! @author H.Tomita and SCALE developpers
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
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  use mod_process, only : &
      PRC_MPIstop
  use gtool_history, only: &
       HistoryInit, &
       HistoryAddVariable, &
       HistoryPutAxis, &
       HistoryPutAdditionalAxis, &
       HistoryPut, &
       HistoryGet
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
  !++ included parameters
  !
  include "inc_precision.h"
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------------
  subroutine HIST_setup
    use mod_stdio, only: &
       IO_FID_CONF

    ! only for register
    call TIME_rapstart('FILE I')
    call TIME_rapend  ('FILE I')
    call TIME_rapstart('FILE O')
    call TIME_rapend  ('FILE O')

    call TIME_rapstart('FILE O')

    call HistoryInit( &
       'SCALE3 HISTORY OUTPUT', 'SCALE-LES ver. 3', 'AICS/RIKEN', &
       (/'x','y','z'/), &
       (/IMAX,JMAX,KMAX/), &
       (/'X', 'Y', 'Z'/), &
       (/'meter','meter','meter'/), &
       (/'REAL4','REAL4','REAL4'/), &
       namelist_fid = IO_FID_CONF &
       )

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_setup

  !-----------------------------------------------------------------------------
  subroutine HIST_reg( &
       itemid, & ! (out)
       item,   & ! (in)
       desc,   & ! (in)
       units,  & ! (in)
       ndim,   & ! (in)
       xdim,   & ! (in)
       ydim,   & ! (in)
       zdim    & ! (in)
       )
    use mod_process, only: &
         PRC_master, &
         PRC_myrank
    implicit none

    integer,          intent(out) :: itemid
    character(len=*), intent( in) :: item
    character(len=*), intent( in) :: desc
    character(len=*), intent( in) :: units
    integer,          intent( in) :: ndim
    character(len=*), intent( in), optional :: xdim
    character(len=*), intent( in), optional :: ydim
    character(len=*), intent( in), optional :: zdim

    logical :: existed

    character(len=2) :: dims(3)

    dims(1) = 'x'
    if ( present(xdim) ) then
       if ( xdim=='half' ) dims(1) = 'xh'
    end if
    dims(2) = 'y'
    if ( present(ydim) ) then
       if ( ydim=='half' ) dims(2) = 'yh'
    end if
    dims(3) = 'z'
    if ( present(zdim) ) then
       if ( zdim=='half' ) dims(3) = 'zh'
    end if

    call TIME_rapstart('FILE O')

    call HistoryAddVariable(item, dims(1:ndim), desc, units, PRC_master, PRC_myrank, & ! (in)
         itemid = itemid, existed = existed) ! (out)

    if ( .not. existed ) then
       call HIST_put_axes
    end if


    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_reg

  !-----------------------------------------------------------------------------
  subroutine HIST_put_1D( &
      itemid, & ! (in)
      var,    & ! (in)
      dt      & ! (in)
      )
    implicit none

    integer,  intent(in) :: itemid
    real(RP), intent(in) :: var(:)
    real(DP), intent(in) :: dt

    real(RP) :: var2(KMAX)
    integer :: k

    call TIME_rapstart('FILE O')

    do k = 1, KMAX
       var2(k) = var(KS+k-1)
    end do
    call HistoryPut(itemid, var2, dt)

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_put_1D
  !-----------------------------------------------------------------------------
  subroutine HIST_put_2D( &
      itemid, & ! (in)
      var,    & ! (in)
      dt      & ! (in)
      )
    implicit none

    integer,  intent(in) :: itemid
    real(RP), intent(in) :: var(:,:)
    real(DP), intent(in) :: dt

    real(RP) :: var2(IMAX*JMAX)
    integer :: i, j

    call TIME_rapstart('FILE O')

    do j = 1, JMAX
    do i = 1, IMAX
       var2(i + (j-1)*IMAX) = var(IS+i-1,JS+j-1)
    end do
    end do

    call HistoryPut(itemid, var2, dt)

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_put_2D
  !-----------------------------------------------------------------------------
  subroutine HIST_put_3D( &
      itemid, & ! (in)
      var,    & ! (in)
      dt      & ! (in)
      )
    implicit none

    integer,  intent(in) :: itemid
    real(RP), intent(in) :: var(:,:,:)
    real(DP), intent(in) :: dt

    real(RP) :: var2(IMAX*JMAX*KMAX)
    integer :: i, j, k
    integer :: s(3)
    intrinsic shape

    call TIME_rapstart('FILE O')

    s = shape(var)
    if ( s(1) == 1 ) then
       do j = 1, JMAX
          do i = 1, IMAX
             var2(i + (j-1)*IMAX) = var(1,IS+i-1,JS+j-1)
          end do
       end do
       call HistoryPut(itemid, var2(1:IMAX*JMAX), dt)
    else
       do k = 1, KMAX
          do j = 1, JMAX
             do i = 1, IMAX
                var2(i + (j-1)*IMAX + (k-1)*JMAX*IMAX) = var(KS+k-1,IS+i-1,JS+j-1)
             end do
          end do
       end do
       call HistoryPut(itemid, var2, dt)
    end if

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_put_3D

  !-----------------------------------------------------------------------------
  subroutine HIST_in_1D( &
       var,    & ! (in)
       item,   & ! (in)
       desc,   & ! (in)
       units,  & ! (in)
       dt,     & ! (in)
       zdim    & ! (in)
       )
    implicit none

    real(RP),         intent(in) :: var(:)
    character(len=*), intent(in) :: item
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: units
    real(DP),         intent(in) :: dt
    character(len=*), intent(in), optional :: zdim

    integer :: itemid

    call TIME_rapstart('FILE O')

    call HIST_reg( &
         itemid,               & ! (out)
         item, desc, units, 1, & ! (in)
         zdim = zdim           & ! (in)
         )

    call HIST_put( itemid, var, dt ) ! (in)

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_in_1D
  !-----------------------------------------------------------------------------
  subroutine HIST_in_2D( &
       var,    & ! (in)
       item,   & ! (in)
       desc,   & ! (in)
       units,  & ! (in)
       dt,     & ! (in)
       xdim,   & ! (in)
       ydim    & ! (in)
       )
    implicit none

    real(RP),         intent(in) :: var(:,:)
    character(len=*), intent(in) :: item
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: units
    real(DP),         intent(in) :: dt
    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim

    integer :: itemid

    call TIME_rapstart('FILE O')

    call HIST_reg( &
         itemid,                  & ! (out)
         item, desc, units, 2,    & ! (in)
         xdim = xdim, ydim = ydim & ! (in)
         )

    call HIST_put( itemid, var, dt ) ! (in)

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_in_2D
  !-----------------------------------------------------------------------------
  subroutine HIST_in_3D( &
       var,    & ! (in)
       item,   & ! (in)
       desc,   & ! (in)
       units,  & ! (in)
       dt,     & ! (in)
       xdim,   & ! (in)
       ydim,   & ! (in)
       zdim    & ! (in)
       )
    implicit none

    real(RP),         intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: item
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: units
    real(DP),         intent(in) :: dt
    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim

    integer :: itemid

    call TIME_rapstart('FILE O')

    call HIST_reg( &
         itemid,                         & ! (out)
         item, desc, units, 3,           & ! (in)
         xdim=xdim, ydim=ydim, zdim=zdim & ! (in)
         )

    call HIST_put( itemid, var, dt ) ! (in)

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_in_3D

  !-----------------------------------------------------------------------------
  ! interface HIST_get
  !-----------------------------------------------------------------------------
  subroutine HIST_get_1D( &
       var,       &
       basename,  &
       varname,   &
       step       &
       )
    use mod_process, only: &
       PRC_myrank
    implicit none

    real(RP),         intent(out) :: var(:)
    character(len=*), intent( in) :: basename
    character(len=*), intent( in) :: varname
    integer,          intent( in) :: step

    call TIME_rapstart('FILE I')

    call HistoryGet( var,                                   & ! (out)
         basename, varname, step, PRC_myrank, single=.true. ) ! (in)

    call TIME_rapend  ('FILE I')

    return
  end subroutine HIST_get_1D
  subroutine HIST_get_2D( &
       var,          &
       basename,     &
       varname,      &
       step,         &
       allow_missing &
       )
    use mod_process, only: &
       PRC_myrank
    implicit none

    real(RP),         intent(out) :: var(:,:)
    character(len=*), intent( in) :: basename
    character(len=*), intent( in) :: varname
    integer,          intent( in) :: step
    logical,          intent( in), optional :: allow_missing

    call TIME_rapstart('FILE I')

    call HistoryGet( var,                                   & ! (out)
         basename, varname, step, PRC_myrank, allow_missing ) ! (in)

    call TIME_rapend  ('FILE I')

    return
  end subroutine HIST_get_2D
  subroutine HIST_get_3D( &
       var,          &
       basename,     &
       varname,      &
       step,         &
       allow_missing &
       )
    use mod_process, only: &
       PRC_myrank
    implicit none

    real(RP),         intent(out) :: var(:,:,:)
    character(len=*), intent( in) :: basename
    character(len=*), intent( in) :: varname
    integer,          intent( in) :: step
    logical,          intent( in), optional :: allow_missing

    call TIME_rapstart('FILE I')

    call HistoryGet( var,                                   & ! (out)
         basename, varname, step, PRC_myrank, allow_missing ) ! (in)

    call TIME_rapend  ('FILE I')

    return
  end subroutine HIST_get_3D

  !-----------------------------------------------------------------------------
  subroutine HIST_write
    use gtool_history, only : &
         HistoryWriteAll
    use mod_time, only : &
         TIME_NOWSEC
    implicit none

    call TIME_rapstart('FILE O')

    call HistoryWriteAll( TIME_NOWSEC )

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_write

  !-----------------------------------------------------------------------------
  ! private
  !-----------------------------------------------------------------------------
  subroutine HIST_put_axes
    use mod_grid, only : &
         GRID_CZ, &
         GRID_CX, &
         GRID_CY, &
         GRID_FZ, &
         GRID_FX, &
         GRID_FY, &
         GRID_CDZ, &
         GRID_CDX, &
         GRID_CDY, &
         GRID_FDZ, &
         GRID_FDX, &
         GRID_FDY, &
         GRID_CBFZ, &
         GRID_CBFX, &
         GRID_CBFY, &
         GRID_FBFZ, &
         GRID_FBFX, &
         GRID_FBFY
    use gtool_history, only : &
         HistoryPutAxis
    implicit none

    call TIME_rapstart('FILE O')

    call HistoryPutAxis('x', GRID_CX(IS:IE))
    call HistoryPutAxis('y', GRID_CY(JS:JE))
    call HistoryPutAxis('z', GRID_CZ(KS:KE))

    call HistoryPutAdditionalAxis('xh', 'X (half level)', 'm', 'xh', GRID_FX(IS:IE))
    call HistoryPutAdditionalAxis('yh', 'Y (half level)', 'm', 'yh', GRID_FY(JS:JE))
    call HistoryPutAdditionalAxis('zh', 'Z (half level)', 'm', 'zh', GRID_FZ(KS:KE))

    call HistoryPutAdditionalAxis('CZ', 'Grid Center Position Z', 'm', 'CZ', GRID_CZ)
    call HistoryPutAdditionalAxis('CX', 'Grid Center Position X', 'm', 'CX', GRID_CX)
    call HistoryPutAdditionalAxis('CY', 'Grid Center Position Y', 'm', 'CY', GRID_CY)
    call HistoryPutAdditionalAxis('FZ', 'Grid Face Position Z', 'm', 'FZ', GRID_FZ)
    call HistoryPutAdditionalAxis('FX', 'Grid Face Position X', 'm', 'FX', GRID_FX)
    call HistoryPutAdditionalAxis('FY', 'Grid Face Position Y', 'm', 'FY', GRID_FY)

    call HistoryPutAdditionalAxis('CDZ', 'Grid Cell length Z', 'm', 'CZ', GRID_CDZ)
    call HistoryPutAdditionalAxis('CDX', 'Grid Cell length X', 'm', 'CX', GRID_CDX)
    call HistoryPutAdditionalAxis('CDY', 'Grid Cell length Y', 'm', 'CY', GRID_CDY)
    call HistoryPutAdditionalAxis('FDZ', 'Grid distance Z', 'm', 'FDZ', GRID_FDZ)
    call HistoryPutAdditionalAxis('FDX', 'Grid distance X', 'm', 'FDX', GRID_FDX)
    call HistoryPutAdditionalAxis('FDY', 'Grid distance Y', 'm', 'FDY', GRID_FDY)

    call HistoryPutAdditionalAxis('CBFZ', 'Boundary factor Center Z', '1', 'CZ', GRID_CBFZ)
    call HistoryPutAdditionalAxis('CBFX', 'Boundary factor Center X', '1', 'CX', GRID_CBFX)
    call HistoryPutAdditionalAxis('CBFY', 'Boundary factor Center Y', '1', 'CY', GRID_CBFY)
    call HistoryPutAdditionalAxis('FBFZ', 'Boundary factor Face Z', '1', 'CZ', GRID_FBFZ)
    call HistoryPutAdditionalAxis('FBFX', 'Boundary factor Face X', '1', 'CX', GRID_FBFX)
    call HistoryPutAdditionalAxis('FBFY', 'Boundary factor Face Y', '1', 'CY', GRID_FBFY)

    call TIME_rapend  ('FILE O')

    return
  end subroutine HIST_put_axes

end module mod_history
!-------------------------------------------------------------------------------
