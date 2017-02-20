!-------------------------------------------------------------------------------------------
!> module NET2G calender
!!
!! @par Description
!!          Calender module for post-process of scale
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_net2g_calender
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_net2g_vars
  use mod_net2g_error

  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: cal_init
  public :: cal_increment

  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: cal_inc_sec
  private :: cal_inc_min
  private :: cal_inc_hour
  private :: cal_date

  !-----------------------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------------------

  !> calender initialization
  !-----------------------------------------------------------------------------------------
  subroutine cal_init( &
      yy,              & ! [inout]
      mm,              & ! [inout]
      dd,              & ! [inout]
      hh,              & ! [inout]
      mn,              & ! [inout]
      sc,              & ! [inout]
      DELT             ) ! [out  ]
    implicit none

    integer, intent(inout) :: yy, mm, dd
    integer, intent(inout) :: hh, mn, sc
    character(*), intent(out) :: DELT

    real(DP) :: inc
    integer  :: i, tint
    character(2) :: tunit
    character(3) :: cint
    !---------------------------------------------------------------------------

    yy = TIME_STARTDATE(1)
    mm = TIME_STARTDATE(2)
    dd = TIME_STARTDATE(3)
    hh = TIME_STARTDATE(4)
    mn = TIME_STARTDATE(5)
    sc = TIME_STARTDATE(6)

    if ( START_TSTEP > 1 ) then
       do i=1, START_TSTEP-1
          call cal_increment( yy, mm, dd, hh, mn, sc )
       enddo
    endif

    inc = HISTORY_DEFAULT_TINTERVAL * dble(INC_TSTEP)
    select case ( HISTORY_DEFAULT_TUNIT )
    case ( "SEC", "sec" )
       if ( inc < 60.0D0 ) then
          if ( LOUT ) write( FID_LOG, '(1X,A)') &
                      "*** WARNING: HISTORY_DEFAULT_TINTERVAL is not compatible!"
          if ( LOUT ) write( FID_LOG, '(1X,A,I7,A)') &
                      "*** ", int(HISTORY_DEFAULT_TINTERVAL), " is too short for Grads"
          tint  = 1     !tentative
          tunit = "mn"  !tentative
       else
          inc = inc / 60.0D0
          if ( inc < 60.0D0 ) then
             tint  = int(inc)
             tunit = "mn"
          else
             inc = inc / 60.0D0
             tint  = int(inc)
             tunit = "hr"
          endif
       endif
    case ( "MIN", "min" )
       if ( inc < 60.0D0 ) then
          tint  = int(inc)
          tunit = "mn"
       else
          inc = inc / 60.0D0
          tint  = int(inc)
          tunit = "hr"
       endif
    case ( "HOUR", "hour" )
       tint  = int(inc)
       tunit = "hr"
    case default
        call err_abort( 0, __LINE__, loc_cal )
    end select

    if ( tint <= 0 ) then
       tint = 1  ! avoid zero
    endif
    write(cint,'(I3)') tint
    DELT = trim(cint)//tunit

    return
  end subroutine cal_init


  !> calender increment calculation
  !-----------------------------------------------------------------------------------------
  subroutine cal_increment( &
      yy,        & ! [inout]
      mm,        & ! [inout]
      dd,        & ! [inout]
      hh,        & ! [inout]
      mn,        & ! [inout]
      sc         ) ! [inout]
    implicit none

    integer, intent(inout) :: yy, mm, dd
    integer, intent(inout) :: hh, mn, sc

    real(DP) :: inc
    !---------------------------------------------------------------------------

    inc = HISTORY_DEFAULT_TINTERVAL * dble(INC_TSTEP)

    select case ( HISTORY_DEFAULT_TUNIT )
    case ( "SEC", "sec" )
       if ( inc < 60.0D0 ) then
          call cal_inc_sec( int(inc), yy, mm, dd, hh, mn, sc )
       else
          inc = inc / 60.0D0
          if ( inc < 60.0D0 ) then
             call cal_inc_min( int(inc), yy, mm, dd, hh, mn )
          else
             inc = inc / 60.0D0
             call cal_inc_hour( int(inc), yy, mm, dd, hh )
          endif
       endif
    case ( "MIN", "min" )
       if ( inc < 60.0D0 ) then
          call cal_inc_min( int(inc), yy, mm, dd, hh, mn )
       else
          inc = inc / 60.0D0
          call cal_inc_hour( int(inc), yy, mm, dd, hh )
       endif
    case ( "HOUR", "hour" )
       call cal_inc_hour( int(inc), yy, mm, dd, hh )
    case default
        call err_abort( 0, __LINE__, loc_cal )
    end select

    return
  end subroutine cal_increment


  !> calender increment calculation: for [sec] increment
  !-----------------------------------------------------------------------------------------
  subroutine cal_inc_sec( &
      inc,    & ! [in]
      yy,     & ! [inout]
      mm,     & ! [inout]
      dd,     & ! [inout]
      hh,     & ! [inout]
      mn,     & ! [inout]
      sc      ) ! [inout]
    implicit none

    integer, intent(in)    :: inc
    integer, intent(inout) :: yy, mm, dd
    integer, intent(inout) :: hh, mn, sc

    integer  :: eday
    !---------------------------------------------------------------------------

    call cal_date( yy, mm, eday )

    sc = sc + inc
    if ( sc >= 60 ) then
       mn = mn + 1
       sc = sc - 60
       if ( mn >= 60 ) then
          hh = hh + 1
          mn = mn - 60
          if ( hh >= 24 ) then
             dd = dd + 1
             hh = hh - 24
             if ( dd > eday ) then
                mm = mm + 1
                dd = dd - eday
                if ( mm > 12 ) then
                   yy = yy + 1
                   mm = mm - 12
                endif
             endif
          endif
       endif
    endif

    return
  end subroutine cal_inc_sec


  !> calender increment calculation: for [min] increment
  !-----------------------------------------------------------------------------------------
  subroutine cal_inc_min( &
      inc,    & ! [in]
      yy,     & ! [inout]
      mm,     & ! [inout]
      dd,     & ! [inout]
      hh,     & ! [inout]
      mn      ) ! [inout]
    implicit none

    integer, intent(in)    :: inc
    integer, intent(inout) :: yy, mm, dd
    integer, intent(inout) :: hh, mn

    integer  :: eday
    !---------------------------------------------------------------------------

    call cal_date( yy, mm, eday )

    mn = mn + inc
    if ( mn >= 60 ) then
       hh = hh + 1
       mn = mn - 60
       if ( hh >= 24 ) then
          dd = dd + 1
          hh = hh - 24
          if ( dd > eday ) then
             mm = mm + 1
             dd = dd - eday
             if ( mm > 12 ) then
                yy = yy + 1
                mm = mm - 12
             endif
          endif
       endif
    endif

    return
  end subroutine cal_inc_min


  !> calender increment calculation: for [hour] increment
  !-----------------------------------------------------------------------------------------
  subroutine cal_inc_hour( &
      inc,    & ! [in]
      yy,     & ! [inout]
      mm,     & ! [inout]
      dd,     & ! [inout]
      hh      ) ! [inout]
    implicit none

    integer, intent(in)    :: inc
    integer, intent(inout) :: yy, mm, dd
    integer, intent(inout) :: hh

    integer  :: eday
    !---------------------------------------------------------------------------

    call cal_date( yy, mm, eday )

    hh = hh + inc
    if ( hh >= 24 ) then
       dd = dd + 1
       hh = hh - 24
       if ( dd > eday ) then
          mm = mm + 1
          dd = dd - eday
          if ( mm > 12 ) then
             yy = yy + 1
             mm = mm - 12
          endif
       endif
    endif

    return
  end subroutine cal_inc_hour


  !> calender date calculation
  !-----------------------------------------------------------------------------------------
  subroutine cal_date( &
      yy,  & ! [in ]
      mm,  & ! [in ]
      dd   ) ! [out]
    implicit none

    integer, intent(in)  :: yy
    integer, intent(in)  :: mm
    integer, intent(out) :: dd

    integer :: rem4, rem100, rem400
    !---------------------------------------------------------------------------

    select case ( mm )
    case ( 4, 6, 9, 11 )
       dd = 30
    case ( 1, 3, 5, 7, 8, 10, 12 )
       dd = 31
    case ( 2 )
       rem4   = int( mod(real(yy), 4.0  ) )
       rem100 = int( mod(real(yy), 100.0) )
       rem400 = int( mod(real(yy), 400.0) )
       dd = 28
       if ( rem4 == 0 ) then            ! T -> leap year
          if ( rem100 == 0 ) then       ! F -> leap year
             if ( rem400 == 0 ) dd = 29 ! T -> leap year
          else
             dd = 29
          endif
       endif
    case ( 0 ) ! ideal test case
       dd = 2147483647
    case default
        call err_abort( 0, __LINE__, loc_cal )
    end select

    return
  end subroutine cal_date

end module mod_net2g_calender
