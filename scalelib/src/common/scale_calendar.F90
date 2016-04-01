!-------------------------------------------------------------------------------
!> module CALENDAR
!!
!! @par Description
!!          gregorian calendar module
!!          this module is available in gregorian calendar date
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-01-29 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_calendar
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
  public :: CALENDAR_setup
  public :: CALENDAR_getDayOfYear
  public :: CALENDAR_date2daysec
  public :: CALENDAR_daysec2date
  public :: CALENDAR_ymd2absday
  public :: CALENDAR_hms2abssec
  public :: CALENDAR_adjust_daysec
  public :: CALENDAR_combine_daysec
  public :: CALENDAR_unit2sec
  public :: CALENDAR_CFunits2sec
  public :: CALENDAR_date2char

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
  real(DP), private :: CALENDAR_DOI  = 365.0_DP !< days of year
  real(DP), private :: CALENDAR_HOUR = 24.0_DP  !< hours   of day
  real(DP), private :: CALENDAR_MIN  = 60.0_DP  !< minutes of hour
  real(DP), private :: CALENDAR_SEC  = 60.0_DP  !< seconds of minute

  integer,  private, parameter :: I_nonleapyear = 1   !< [index] non leap year
  integer,  private, parameter :: I_leapyear    = 2   !< [index] leap year
  integer,  private            :: dayofmonth(12,2)    !< days of each month
  data dayofmonth / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, & ! non-leap year
                    31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  / ! leap year

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CALENDAR_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CALENDAR] / Categ[COMMON] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** No namelists.'

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
    !---------------------------------------------------------------------------

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

    gmonth_mod = mod( gmonth-1, 12 ) + 1
    gyear_mod  = gyear + ( gmonth-gmonth_mod ) / 12

    yearday = int( CALENDAR_DOI * gyear_mod )                   &
            + int( real(gyear_mod+oyear-1,kind=DP) /   4.0_DP ) &
            - int( real(gyear_mod+oyear-1,kind=DP) / 100.0_DP ) &
            + int( real(gyear_mod+oyear-1,kind=DP) / 400.0_DP ) &
            - int( real(          oyear-1,kind=DP) /   4.0_DP ) &
            + int( real(          oyear-1,kind=DP) / 100.0_DP ) &
            - int( real(          oyear-1,kind=DP) / 400.0_DP )

    ileap = I_nonleapyear
    if( checkleap(gyear_mod+oyear) ) ileap = I_leapyear

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

    gyear = int( real(absday,kind=DP) / 366.0_DP ) ! first guess
    do i = 1, 1000
       call CALENDAR_ymd2absday( checkday, gyear+1, 1, 1, oyear )
       if( absday < checkday ) exit
       gyear = gyear + 1
    enddo

    ileap = I_nonleapyear
    if( checkleap(gyear+oyear) ) ileap = I_leapyear

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
    use scale_process, only: &
       PRC_MPIstop
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
       write(*,*) ' xxx Unsupported UNIT: ', trim(unit), ', ', value
       call PRC_MPIstop
    endselect

    return
  end subroutine CALENDAR_unit2sec

  !-----------------------------------------------------------------------------
  !> Convert time in units of the CF convention to second
  function CALENDAR_CFunits2sec( cftime, cfunits, offset_year, startdaysec ) result( sec )
    use scale_process, only: &
       PRC_MPIstop
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

    l = index(cfunits, " since ")
    if ( l > 1 ) then ! untis is under the CF convension
       tunit = cfunits(1:l-1)
       buf = cfunits(l+7:)

       l = index(buf, "-")
       if ( l .ne. 5 ) then
          write(*,*) 'xxx units for time is invalid (year): ', trim(cfunits), ' ', trim(buf)
          call PRC_MPIstop
       end if
       read( buf(1:4), *) date(1) ! year
       date(1) = date(1) - offset_year
       buf = buf(6:)

       l = index(buf, "-")
       if ( l .ne. 3 ) then
          write(*,*) 'xxx units for time is invalid (month): ', trim(cfunits), ' ', trim(buf)
          call PRC_MPIstop
       end if
       read( buf(1:2), *) date(2) ! month
       buf = buf(4:)

       l = index(buf, " ")
       if ( l .ne. 3 ) then
          write(*,*) 'xxx units for time is invalid (day): ', trim(cfunits), ' ', trim(buf)
          call PRC_MPIstop
       end if
       read( buf(1:2), *) date(3) ! day
       buf = buf(4:)

       l = index(buf, ":")
       if ( l .ne. 3 ) then
          write(*,*) 'xxx units for time is invalid (hour): ', trim(cfunits), ' ', trim(buf)
          call PRC_MPIstop
       end if
       read( buf(1:2), *) date(4) ! hour
       buf = buf(4:)

       l = index(buf, ":")
       if ( l .ne. 3 ) then
          write(*,*) 'xxx units for time is invalid (min): ', trim(cfunits), ' ', trim(buf)
          call PRC_MPIstop
       end if
       read( buf(1:2), *) date(5) ! min
       buf = buf(4:)

       if ( len_trim(buf) .ne. 2 ) then
          write(*,*) 'xxx units for time is invalid (sec): ', trim(cfunits), ' ', trim(buf), len_trim(buf)
          call PRC_MPIstop
       end if
       read( buf(1:2), *) date(6) ! sec

       call CALENDAR_date2daysec( day,     & ! (out)
                                  sec0,    & ! (out)
                                  date(:), & ! (in)
                                  0.0_DP,  & ! (in)
                                  offset_year ) ! (in)

       sec0 = CALENDAR_combine_daysec( day, sec0 )
    else
       tunit = cfunits
       sec0 = startdaysec
    end if

    call CALENDAR_unit2sec(sec, cftime, tunit)

    sec = sec0 + sec

    return
  end function CALENDAR_CFunits2sec

  !-----------------------------------------------------------------------------
  !> Convert from gregorian date to absolute day/second
  subroutine CALENDAR_date2char( &
       chardate,   &
       ymdhms,     &
       subsec,     &
       offset_year )
    implicit none

    character(len=27), intent(out) :: chardate    !< formatted date character
    integer,           intent(in)  :: ymdhms(6)   !< date
    real(DP),          intent(in)  :: subsec      !< subsecond
    integer,           intent(in)  :: offset_year !< offset year
    !---------------------------------------------------------------------------

    write(chardate,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A2,F6.3)')   &
                       ymdhms(1)+offset_year,'/',ymdhms(2),'/',ymdhms(3),' ',  &
                       ymdhms(4)            ,':',ymdhms(5),':',ymdhms(6),' +', &
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

    mjd_day = int( 365.25_DP * (y+oyear) )                          & ! year
            + int( (y+oyear)/400.0_DP ) - int( (y+oyear)/100.0_DP ) & ! leap year
            + int( 30.59_DP * m-2 )                                 & ! month
            + ymdhms(I_day)                                         & ! day
            + 678912                                                ! constant

    mjd     = real(mjd_day,kind=DP)       & ! day
            + ymdhms(I_hour) /    24.0_DP & ! hour
            + ymdhms(I_min)  /  1440.0_DP & ! min
            + ymdhms(I_sec)  / 86400.0_DP   ! sec

    return
  end subroutine CALENDAR_ymdhms2mjd

end module scale_calendar
!-------------------------------------------------------------------------------
