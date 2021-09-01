!-------------------------------------------------------------------------------
!> module CALENDAR
!!
!! @par Description
!!          gregorian calendar module
!!          this module is available in gregorian calendar date
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_calendar
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CALENDAR_setup
  public :: CALENDAR_getDayOfYear
  public :: CALENDAR_date2daysec
  public :: CALENDAR_daysec2date
  public :: CALENDAR_ymd2absday
  public :: CALENDAR_hms2abssec
  public :: CALENDAR_adjust_daysec
  public :: CALENDAR_combine_daysec
  public :: CALENDAR_unit2sec
  public :: CALENDAR_sec2unit
  public :: CALENDAR_CFunits2sec
  public :: CALENDAR_date2char
  public :: CALENDAR_get_name

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: I_year  = 1 !< [index] year
  integer, public, parameter :: I_month = 2 !< [index] month
  integer, public, parameter :: I_day   = 3 !< [index] day
  integer, public, parameter :: I_hour  = 4 !< [index] hour
  integer, public, parameter :: I_min   = 5 !< [index] minute
  integer, public, parameter :: I_sec   = 6 !< [index] second

  real(DP), public :: CALENDAR_DOI  = 365.0_DP !< days of year
  real(DP), public :: CALENDAR_HOUR = 24.0_DP  !< hours   of day
  real(DP), public :: CALENDAR_MIN  = 60.0_DP  !< minutes of hour
  real(DP), public :: CALENDAR_SEC  = 60.0_DP  !< seconds of minute

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CALENDAR_absday2ymd
  private :: CALENDAR_abssec2hms
  private :: CALENDAR_ymdhms2nd
  private :: CALENDAR_ymdhms2mjd
  private :: checkleap

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private :: CALENDAR_360DAYS = .false.
  logical,  private :: CALENDAR_365DAYS = .false.
  logical,  private :: CALENDAR_USER    = .false.

  integer,  private, parameter :: I_nonleapyear = 1   !< [index] non leap year
  integer,  private, parameter :: I_leapyear    = 2   !< [index] leap year
  integer,  private, parameter :: I_360days     = 3   !< [index] 360 days
  integer,  private, parameter :: I_USER        = 4   !< [index] user defined calendar

  integer,  private            :: dayofmonth(12,4)    !< days of each month
  data dayofmonth / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, & ! non-leap year
                    31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, & ! leap year
                    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, & ! 360 days
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /  ! for CALENDAR_USER_DEFINED

  integer,  private :: CALENDAR_USER_DEFINED(12)
  integer,  private :: CALENDAR_MONTH = 12
  logical,  private :: debug = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CALENDAR_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_CALENDAR / &
       CALENDAR_360DAYS,        &
       CALENDAR_365DAYS,        &
       CALENDAR_HOUR,           &
       CALENDAR_MIN,            &
       CALENDAR_SEC,            &
       CALENDAR_USER_DEFINED,   &
       debug

    integer :: ierr, i
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CALENDAR_setup",*) 'Setup'

    CALENDAR_USER_DEFINED(:) = 0

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CALENDAR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CALENDAR_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CALENDAR_setup",*) 'Not appropriate names in namelist PARAM_CALENDAR. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CALENDAR)

    LOG_NEWLINE
    LOG_INFO("CALENDAR_setup",*) 'Calendar settings'
    if ( maxval(CALENDAR_USER_DEFINED) > 0 ) then
       dayofmonth(:,I_USER) = CALENDAR_USER_DEFINED(:)
       CALENDAR_DOI  = sum(CALENDAR_USER_DEFINED)
       CALENDAR_USER = .true.
       LOG_INFO_CONT(*) 'DayOfYear = ', int(CALENDAR_DOI), ' : user defined calendar'
       do i = 1, 12
          if ( CALENDAR_USER_DEFINED(i) > 0 ) then
             LOG_INFO_CONT(*) 'month #', i, ': ', CALENDAR_USER_DEFINED(i), "days"
          else
             CALENDAR_MONTH = i - 1
             exit
          endif
       enddo
    elseif( CALENDAR_360DAYS ) then
       CALENDAR_DOI = 360.0_DP
       LOG_INFO_CONT(*) 'DayOfYear = 360 : ideal setting'
    elseif( CALENDAR_365DAYS ) then
       CALENDAR_DOI = 365.0_DP
       LOG_INFO_CONT(*) 'DayOfYear = 365 : ideal setting'
    else
       LOG_INFO_CONT(*) 'DayOfYear = 365 or 366 : Gregorian calendar'
    endif

    if (int(CALENDAR_SEC) /= CALENDAR_SEC .or. int(CALENDAR_SEC) == 0 ) then
       LOG_ERROR("CALENDAR_setup",*) 'CALENDAR_SEC must be a natural number!'
       call abort
    elseif (int(CALENDAR_MIN) /= CALENDAR_MIN .or. int(CALENDAR_MIN) == 0 ) then
       LOG_ERROR("CALENDAR_setup",*) 'CALENDAR_MIN must be a natural number!'
       call abort
    elseif (int(CALENDAR_HOUR) /= CALENDAR_HOUR .or. int(CALENDAR_HOUR) == 0 ) then
       LOG_ERROR("CALENDAR_setup",*) 'CALENDAR_HOUR must be a natural number!'
       call abort
    endif


    return
  end subroutine CALENDAR_setup

  !-----------------------------------------------------------------------------
  !> Get day of year
  subroutine CALENDAR_getDayOfYear( &
       DayOfYear, &
       iyear      )
    implicit none

    real(DP), intent(out) :: DayOfYear ! # of day of year
    integer,  intent(in)  :: iyear     ! current year
    !---------------------------------------------------------------------------

    DayOfYear = CALENDAR_DOI
    if( checkleap(iyear) ) DayOfYear = CALENDAR_DOI + 1.0_DP

    return
  end subroutine CALENDAR_getDayOfYear

  !-----------------------------------------------------------------------------
  !> Convert from gregorian date to absolute day/second
  subroutine CALENDAR_date2daysec( &
       absday,     &
       abssec,     &
       ymdhms,     &
       subsec,     &
       offset_year )
    implicit none

    integer,  intent(out) :: absday      !< absolute day
    real(DP), intent(out) :: abssec      !< absolute second
    integer,  intent(in)  :: ymdhms(6)   !< date
    real(DP), intent(in)  :: subsec      !< subsecond
    integer,  intent(in)  :: offset_year !< offset year
    integer               :: ileap
    logical               :: date_error = .false.
    !---------------------------------------------------------------------------

    ! check date
    if ( ymdhms(I_month) < 1 .or. CALENDAR_MONTH < ymdhms(I_month) ) then
       LOG_ERROR("CALENDAR_date2daysec",*) 'Inputted month does not match to the calendar.'
       date_error = .true.
    endif

    ileap = I_nonleapyear
    if( checkleap(ymdhms(I_year)+offset_year) ) ileap = I_leapyear
    if( CALENDAR_360DAYS     ) ileap = I_360days
    if( CALENDAR_USER        ) ileap = I_USER

    if ( ymdhms(I_day ) < 1 .or. dayofmonth(ymdhms(I_month),ileap) < ymdhms(I_day) ) then
       LOG_ERROR("CALENDAR_date2daysec",*) 'Inputted day does not match to the calendar.'
       date_error = .true.
    endif
    if ( ymdhms(I_hour) < 0 .or. (int(CALENDAR_HOUR)-1) < ymdhms(I_hour) ) then
       LOG_ERROR("CALENDAR_date2daysec",*) 'Inputted hour does not match to the calendar.'
       date_error = .true.
    endif
    if ( ymdhms(I_min ) < 0 .or. (int(CALENDAR_MIN )-1) < ymdhms(I_min ) ) then
       LOG_ERROR("CALENDAR_date2daysec",*) 'Inputted minute does not match to the calendar.'
       date_error = .true.
    endif
    if ( ymdhms(I_sec ) < 0 .or. (int(CALENDAR_SEC )-1) < ymdhms(I_sec ) ) then
       LOG_ERROR("CALENDAR_date2daysec",*) 'Inputted second does not match to the calendar.'
       date_error = .true.
    endif
    if( date_error ) call abort

    call CALENDAR_ymd2absday( absday,          & ! [OUT]
                              ymdhms(I_year),  & ! [IN]
                              ymdhms(I_month), & ! [IN]
                              ymdhms(I_day),   & ! [IN]
                              offset_year      ) ! [IN]

    call CALENDAR_hms2abssec( abssec,          & ! [OUT]
                              ymdhms(I_hour),  & ! [IN]
                              ymdhms(I_min),   & ! [IN]
                              ymdhms(I_sec),   & ! [IN]
                              subsec           ) ! [IN]

    return
  end subroutine CALENDAR_date2daysec

  !-----------------------------------------------------------------------------
  !> Convert from gregorian date to absolute day/second
  subroutine CALENDAR_daysec2date( &
       ymdhms,     &
       subsec,     &
       absday,     &
       abssec,     &
       offset_year )
    implicit none

    integer,  intent(out) :: ymdhms(6)   !< date
    real(DP), intent(out) :: subsec      !< subsecond
    integer,  intent(in)  :: absday      !< absolute day
    real(DP), intent(in)  :: abssec      !< absolute second
    integer,  intent(in)  :: offset_year !< offset year
    !---------------------------------------------------------------------------

    call CALENDAR_absday2ymd( ymdhms(I_year),  & ! [OUT]
                              ymdhms(I_month), & ! [OUT]
                              ymdhms(I_day),   & ! [OUT]
                              absday,          & ! [IN]
                              offset_year      ) ! [IN]

    call CALENDAR_abssec2hms( ymdhms(I_hour),  & ! [OUT]
                              ymdhms(I_min),   & ! [OUT]
                              ymdhms(I_sec),   & ! [OUT]
                              subsec,          & ! [OUT]
                              abssec           ) ! [IN]

    return
  end subroutine CALENDAR_daysec2date

  !-----------------------------------------------------------------------------
  !> Convert from gregorian date to absolute day, DAY 0 is AD1/1/1
  subroutine CALENDAR_ymd2absday( &
       absday, &
       gyear,  &
       gmonth, &
       gday,   &
       oyear   )
    implicit none

    integer, intent(out) :: absday !< absolute day
    integer, intent(in)  :: gyear  !< year
    integer, intent(in)  :: gmonth !< month
    integer, intent(in)  :: gday   !< day
    integer, intent(in)  :: oyear  !< offset year

    integer :: gyear_mod, gmonth_mod
    integer :: yearday, monthday
    integer :: m, ileap
    !---------------------------------------------------------------------------

    gmonth_mod = mod( gmonth-1, CALENDAR_MONTH ) + 1
    gyear_mod  = gyear + ( gmonth-gmonth_mod ) / CALENDAR_MONTH

    if ( CALENDAR_360DAYS .OR. CALENDAR_365DAYS .or. CALENDAR_USER) then
       yearday = int( CALENDAR_DOI * ( gyear_mod - oyear ) )
    else
       yearday = int( CALENDAR_DOI * ( gyear_mod - oyear ) ) &
               + int( real(gyear_mod-1,kind=DP) /   4.0_DP ) &
               - int( real(gyear_mod-1,kind=DP) / 100.0_DP ) &
               + int( real(gyear_mod-1,kind=DP) / 400.0_DP ) &
               - int( real(oyear    -1,kind=DP) /   4.0_DP ) &
               + int( real(oyear    -1,kind=DP) / 100.0_DP ) &
               - int( real(oyear    -1,kind=DP) / 400.0_DP )
    endif

    ileap = I_nonleapyear
    if( checkleap(gyear_mod) ) ileap = I_leapyear
    if( CALENDAR_360DAYS     ) ileap = I_360days
    if( CALENDAR_USER        ) ileap = I_USER

    monthday = 0
    do m = 1, gmonth_mod-1
       monthday = monthday + dayofmonth(m,ileap)
    enddo

    absday = yearday + monthday + gday - 1

    return
  end subroutine CALENDAR_ymd2absday

  !-----------------------------------------------------------------------------
  !> Convert from absolute day to gregorian date, DAY 0 is AD1/1/1
  subroutine CALENDAR_absday2ymd( &
       gyear,  &
       gmonth, &
       gday,   &
       absday, &
       oyear   )
    implicit none

    integer, intent(out) :: gyear  !< year
    integer, intent(out) :: gmonth !< month
    integer, intent(out) :: gday   !< day
    integer, intent(in)  :: absday !< absolute day
    integer, intent(in)  :: oyear  !< offset year

    integer :: checkday
    integer :: i, ileap
    !---------------------------------------------------------------------------

    gyear = int( real(absday,kind=DP) / (CALENDAR_DOI+1.0_DP) ) + oyear ! first guess

    do i = 1, 1000
       call CALENDAR_ymd2absday( checkday, gyear+1, 1, 1, oyear )
       if( absday < checkday ) exit
       gyear = gyear + 1
    enddo

    ileap = I_nonleapyear
    if( checkleap(gyear) ) ileap = I_leapyear
    if( CALENDAR_360DAYS ) ileap = I_360days
    if( CALENDAR_USER    ) ileap = I_USER

    gmonth = 1
    do i = 1, 1000
       call CALENDAR_ymd2absday( checkday, gyear, gmonth, dayofmonth(gmonth,ileap), oyear )
       if( absday <= checkday ) exit
       gmonth = gmonth + 1
    enddo

    call CALENDAR_ymd2absday( checkday, gyear, gmonth, 1, oyear )
    gday = absday - checkday + 1

    return
  end subroutine CALENDAR_absday2ymd

  !-----------------------------------------------------------------------------
  !> Hour, minute, second, subsecond -> absolute second
  subroutine CALENDAR_hms2abssec( &
       abssec, &
       hour,   &
       minute, &
       second, &
       subsec  )
    implicit none

    real(DP), intent(out) :: abssec !< absolute second
    integer,  intent(in)  :: hour   !< hour
    integer,  intent(in)  :: minute !< minute
    integer,  intent(in)  :: second !< second
    real(DP), intent(in)  :: subsec !< subsecond
    !---------------------------------------------------------------------------

    abssec = real(hour,  kind=DP) * CALENDAR_MIN * CALENDAR_SEC &
           + real(minute,kind=DP) * CALENDAR_SEC                &
           + real(second,kind=DP)                               &
           + subsec

    return
  end subroutine CALENDAR_hms2abssec

  !-----------------------------------------------------------------------------
  !> Absolute second -> hour, minute, second, subsecond
  subroutine CALENDAR_abssec2hms( &
       hour,   &
       minute, &
       second, &
       subsec, &
       abssec  )
    implicit none

    integer,  intent(out) :: hour   !< hour
    integer,  intent(out) :: minute !< minute
    integer,  intent(out) :: second !< second
    real(DP), intent(out) :: subsec !< subsecond
    real(DP), intent(in)  :: abssec !< absolute second

    real(DP) :: nsec, nmin, nhour, temp
    !---------------------------------------------------------------------------

    nsec   = real(int(abssec),kind=DP)
    subsec = abssec - nsec

    temp   = mod( nsec, CALENDAR_SEC )
    second = int( temp )
    nmin   = ( nsec-temp ) / CALENDAR_SEC

    temp   = mod( nmin, CALENDAR_MIN )
    minute = int( temp )
    nhour  = ( nmin-temp ) / CALENDAR_MIN

    temp   = mod( nhour, CALENDAR_HOUR )
    hour   = int( temp )

    return
  end subroutine CALENDAR_abssec2hms

  !-----------------------------------------------------------------------------
  !> Adjust day and second
  subroutine CALENDAR_adjust_daysec( &
       absday, &
       abssec  )
    implicit none

    integer,  intent(inout) :: absday !< absolute day
    real(DP), intent(inout) :: abssec !< absolute second

    integer :: addday
    !---------------------------------------------------------------------------

    addday = int( abssec / ( CALENDAR_HOUR * CALENDAR_MIN * CALENDAR_SEC ) )

    absday = absday + addday

    abssec = abssec - real(addday,kind=DP) * CALENDAR_HOUR * CALENDAR_MIN * CALENDAR_SEC

    if ( abssec < 0.0_DP ) then
       absday = absday - 1
       abssec = abssec + CALENDAR_HOUR * CALENDAR_MIN * CALENDAR_SEC
    endif

    return
  end subroutine CALENDAR_adjust_daysec

  !-----------------------------------------------------------------------------
  !> Combine day and second
  function CALENDAR_combine_daysec( absday, abssec ) result(daysec)
    implicit none

    integer,  intent(in) :: absday !< absolute day
    real(DP), intent(in) :: abssec !< absolute second
    real(DP)             :: daysec !< absolute day.second
    !---------------------------------------------------------------------------

    daysec = real(absday,kind=DP) * CALENDAR_SEC * CALENDAR_MIN * CALENDAR_HOUR &
           + abssec

    return
  end function CALENDAR_combine_daysec

  !-----------------------------------------------------------------------------
  !> Convert several units to second
  subroutine CALENDAR_unit2sec( &
       second, &
       value,  &
       unit    )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(DP),         intent(out) :: second !< second
    real(DP),         intent(in)  :: value  !< value
    character(len=*), intent(in)  :: unit   !< variable unit
    !---------------------------------------------------------------------------

    select case(unit)
    case('MSEC')
       second = value * 1.E-3_DP
    case('SEC', 'seconds')
       second = value
    case('MIN')
       second = value * CALENDAR_SEC
    case('HOUR')
       second = value * CALENDAR_SEC * CALENDAR_MIN
    case('DAY')
       second = value * CALENDAR_SEC * CALENDAR_MIN * CALENDAR_HOUR
    case default
       LOG_ERROR("CALENDAR_unit2sec",*) 'Unsupported UNIT: ', trim(unit), ', ', value
       call PRC_abort
    endselect

    return
  end subroutine CALENDAR_unit2sec

  !-----------------------------------------------------------------------------
  !> Convert several second to specified unit
  subroutine CALENDAR_sec2unit( &
     value,  &
     second, &
     unit    )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(DP),         intent(out) :: value
    real(DP),         intent( in) :: second
    character(len=*), intent( in) :: unit
    !---------------------------------------------------------------------------

    select case(trim(unit))
    case('MSEC', 'msec')
       value = second / 1.0E-3_DP
    case('SEC', 'seconds', 'sec', 's')
       value = second
    case('MIN', 'mins', 'min')
       value = second / CALENDAR_SEC
    case('HOUR', 'hours', 'hour', 'h')
       value = second / (CALENDAR_SEC * CALENDAR_MIN)
    case('DAY', 'days', 'day')
       value = second / (CALENDAR_SEC * CALENDAR_MIN * CALENDAR_HOUR)
    case default
       LOG_ERROR("CALENDAR_sec2unit",*) 'Unsupported UNIT: ', trim(unit), ', ', value
       call PRC_abort
    endselect

  end subroutine CALENDAR_sec2unit

  !-----------------------------------------------------------------------------
  !> Convert time in units of the CF convention to second
  function CALENDAR_CFunits2sec( cftime, cfunits, offset_year, startdaysec ) result( sec )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(DP),         intent(in) :: cftime
    character(len=*), intent(in) :: cfunits
    integer,          intent(in) :: offset_year
    real(DP),         intent(in), optional :: startdaysec
    real(DP)                     :: sec

    character(len=H_MID) :: tunit
    character(len=H_MID) :: buf

    integer  :: date(6)
    integer  :: day
    real(DP) :: sec0

    integer :: l

    intrinsic index
    !---------------------------------------------------------------------------

    l = index( cfunits, " since " )
    if ( l > 1 ) then ! untis is under the CF convension
       tunit = cfunits(1:l-1)
       buf   = cfunits(l+7:)

       l = index(buf,"-")
       if ( l /= 5 ) then
          LOG_ERROR("CALENDAR_CFunits2sec",*) 'units for time is invalid (year)'
          LOG_ERROR_CONT(*) trim(cfunits)
          LOG_ERROR_CONT(*) trim(buf)
          call PRC_abort
       end if
       read(buf(1:4),*) date(1) ! year
       buf = buf(6:)

       l = index(buf,"-")
       if ( l /= 3 ) then
          LOG_ERROR("CALENDAR_CFunits2sec",*) 'units for time is invalid (month)'
          LOG_ERROR_CONT(*) trim(cfunits)
          LOG_ERROR_CONT(*) trim(buf)
          call PRC_abort
       end if
       read(buf(1:2),*) date(2) ! month
       buf = buf(4:)

       l = index(buf," ")
       if ( l /= 3 ) then
          LOG_ERROR("CALENDAR_CFunits2sec",*) 'units for time is invalid (day)'
          LOG_ERROR_CONT(*) trim(cfunits)
          LOG_ERROR_CONT(*) trim(buf)
          call PRC_abort
       end if
       read(buf(1:2),*) date(3) ! day
       buf = buf(4:)

       l = index(buf,":")
       if ( l /= 3 ) then
          LOG_ERROR("CALENDAR_CFunits2sec",*) 'units for time is invalid (hour)'
          LOG_ERROR_CONT(*) trim(cfunits)
          LOG_ERROR_CONT(*) trim(buf)
          call PRC_abort
       end if
       read(buf(1:2),*) date(4) ! hour
       buf = buf(4:)

       l = index(buf,":")
       if ( l /= 3 ) then
          LOG_ERROR("CALENDAR_CFunits2sec",*) 'units for time is invalid (min)'
          LOG_ERROR_CONT(*) trim(cfunits)
          LOG_ERROR_CONT(*) trim(buf)
          call PRC_abort
       end if
       read(buf(1:2),*) date(5) ! min
       buf = buf(4:)

       if ( len_trim(buf) /= 2 ) then
          LOG_ERROR("CALENDAR_CFunits2sec",*) 'units for time is invalid (sec)'
          LOG_ERROR_CONT(*) trim(cfunits)
          LOG_ERROR_CONT(*) trim(buf)
          LOG_ERROR_CONT(*) len_trim(buf)
          call PRC_abort
       end if
       read(buf(1:2),*) date(6) ! sec

       call CALENDAR_date2daysec( day,        & ! (out)
                                  sec0,       & ! (out)
                                  date(:),    & ! (in)
                                  0.0_DP,     & ! (in)
                                  offset_year ) ! (in)

       sec0 = CALENDAR_combine_daysec( day, sec0 )
    else
       tunit = cfunits
       if ( present(startdaysec) ) then
          sec0 = startdaysec
       else
          sec0 = 0.0_DP
       end if
    end if

    call CALENDAR_unit2sec( sec, cftime, tunit )

    sec = sec0 + sec

    return
  end function CALENDAR_CFunits2sec

  !-----------------------------------------------------------------------------
  !> Convert from gregorian date to absolute day/second
  subroutine CALENDAR_date2char( &
       chardate, &
       ymdhms,   &
       subsec    )
    implicit none

    character(len=27), intent(out) :: chardate    !< formatted date character
    integer,           intent(in)  :: ymdhms(6)   !< date
    real(DP),          intent(in)  :: subsec      !< subsecond

    character(len=4)               :: seclen      !< length of seconds
    character(len=4)               :: minlen      !< length of minutes
    character(len=4)               :: hourlen     !< length of hours
    character(len=4)               :: daylen      !< length of days
    character(len=4)               :: monthlen    !< length of months


    !---------------------------------------------------------------------------
    if ( CALENDAR_USER) then
       write(daylen, '(A1,I1.1,A1,I1.1)') 'I', int(log10(real(maxval(CALENDAR_USER_DEFINED))))+1, &
                                          '.', int(log10(real(maxval(CALENDAR_USER_DEFINED))))+1
       if( CALENDAR_MONTH < 10 ) then
          write(monthlen, '(A4)') 'I1.1'
       else
          write(monthlen, '(A4)') 'I2.2'
       endif
    else
       write(daylen,   '(A4)') 'I2.2'
       write(monthlen, '(A4)') 'I2.2'
    endif

    write(seclen, '(A1,I1.1,A1,I1.1)') 'I', int(log10(max(CALENDAR_SEC -1.0_DP, 1.0_DP)))+1, &
                                       '.', int(log10(max(CALENDAR_SEC -1.0_DP, 1.0_DP)))+1
    write(minlen, '(A1,I1.1,A1,I1.1)') 'I', int(log10(max(CALENDAR_MIN -1.0_DP, 1.0_DP)))+1, &
                                       '.', int(log10(max(CALENDAR_MIN -1.0_DP, 1.0_DP)))+1
    write(hourlen,'(A1,I1.1,A1,I1.1)') 'I', int(log10(max(CALENDAR_HOUR-1.0_DP, 1.0_DP)))+1, &
                                       '.', int(log10(max(CALENDAR_HOUR-1.0_DP, 1.0_DP)))+1

!    write(chardate,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A2,F6.3)')   &
    write(chardate,'(I4.4,A1,'//monthlen//',A1,'//daylen//',A1,'//hourlen//',A1,'//minlen//',A1,'//seclen//',A2,F6.3)') &
                       ymdhms(1),'/',ymdhms(2),'/',ymdhms(3),' ',  &
                       ymdhms(4),':',ymdhms(5),':',ymdhms(6),' +', &
                       subsec

    return
  end subroutine CALENDAR_date2char

  !-----------------------------------------------------------------------------
  !> Check leap year
  !> @return checkleap
  function checkleap( iyear )
    implicit none

    integer, intent(in) :: iyear     !< current year
    logical             :: checkleap

    integer :: check4, check100, check400
    !---------------------------------------------------------------------------

    check4   = mod(iyear,4  )
    check100 = mod(iyear,100)
    check400 = mod(iyear,400)

    checkleap = .false.
    if( check4   == 0 ) checkleap = .true.
    if( check100 == 0 ) checkleap = .false.
    if( check400 == 0 ) checkleap = .true.

    if( CALENDAR_360DAYS ) checkleap = .false.
    if( CALENDAR_365DAYS ) checkleap = .false.
    if( CALENDAR_USER    ) checkleap = .false.

  end function checkleap

  !-----------------------------------------------------------------------------
  !> Calc number of day
  subroutine CALENDAR_ymdhms2nd( &
       nd,     &
       ymdhms, &
       oyear   )
    implicit none

    real(DP), intent(out) :: nd        !< # of day from Jan 1st
    integer,  intent(in)  :: ymdhms(6) !< date
    integer,  intent(in)  :: oyear     !< offset year

    integer :: absday, absday_jan1
    !---------------------------------------------------------------------------

    call CALENDAR_ymd2absday( absday,          & ! [OUT]
                              ymdhms(I_year),  & ! [IN]
                              ymdhms(I_month), & ! [IN]
                              ymdhms(I_day),   & ! [IN]
                              oyear            ) ! [IN]

    call CALENDAR_ymd2absday( absday_jan1,     & ! [OUT]
                              ymdhms(I_year),  & ! [IN]
                              1,               & ! [IN]
                              1,               & ! [IN]
                              oyear            ) ! [IN]

    nd = absday - absday_jan1 + 1.0_DP

    return
  end subroutine CALENDAR_ymdhms2nd

  !-----------------------------------------------------------------------------
  !> Calc modified julian day number (MJD), epoch time is 1858/11/17 00:00:00 UT
  subroutine CALENDAR_ymdhms2mjd( &
       mjd,    &
       ymdhms, &
       oyear   )
    implicit none

    real(DP), intent(out) :: mjd       !< modified julian day number (MJD)
    integer,  intent(in)  :: ymdhms(6) !< date in gregorian calendar
    integer,  intent(in)  :: oyear     !< offset year

    integer :: y, m, mjd_day
    !---------------------------------------------------------------------------

    if (      ymdhms(I_month) == 1 &
         .OR. ymdhms(I_month) == 2 ) then
       y = ymdhms(I_year)  - 1
       m = ymdhms(I_month) + 2
    else
       y = ymdhms(I_year)
       m = ymdhms(I_month)
    endif

    mjd_day = int( 365.25_DP * y )                  & ! year
            + int( y/400.0_DP ) - int( y/100.0_DP ) & ! leap year
            + int( 30.59_DP * m-2 )                 & ! month
            + ymdhms(I_day)                         & ! day
            + 678912                                ! constant

    mjd     = real(mjd_day,kind=DP)                                      & ! day
            + ymdhms(I_hour) / (CALENDAR_HOUR                          ) & ! hour
            + ymdhms(I_min)  / (CALENDAR_HOUR*CALENDAR_MIN             ) & ! min
            + ymdhms(I_sec)  / (CALENDAR_HOUR*CALENDAR_MIN*CALENDAR_SEC)   ! sec

    return
  end subroutine CALENDAR_ymdhms2mjd

  subroutine CALENDAR_get_name(name)
    character(len=*), intent(out) :: name

    if    ( CALENDAR_360DAYS ) then
       name = "360_day"
    elseif( CALENDAR_365DAYS ) then
       name = "365_day"
    elseif( CALENDAR_USER    ) then
       name = "USER_DEFINED"
    else
       name = "gregorian"
    endif

    return
  end subroutine CALENDAR_get_name

end module scale_calendar
