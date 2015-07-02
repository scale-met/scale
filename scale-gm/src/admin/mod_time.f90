!-------------------------------------------------------------------------------
!> Module time management
!!
!! @par Description
!!          This module is for the time management
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_time
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TIME_setup
  public :: TIME_report
  public :: TIME_advance

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=ADM_NSYS), public :: TIME_INTEG_TYPE = 'UNDEF'  ! Integration method in large steps
  !                                                  = 'RK2'    ! Runge-Kutta 2nd
  !                                                  = 'RK3'    ! Runge-Kutta 3rd
  !                                                  = 'RK4'    ! Runge-Kutta 4th
  !                                                  = 'TRCADV' ! Tracer advection only

  logical,  public :: TIME_SPLIT     = .true. ! Horizontally splitting?

  integer,  public :: TIME_LSTEP_MAX = 10     ! Max steps of large step
  integer,  public :: TIME_SSTEP_MAX          ! Max steps of small step

  real(RP), public :: TIME_DTL      = 5.0_RP  ! Time interval for large step [sec]
  real(RP), public :: TIME_DTS                ! Time interval for small step [sec]
  !
  real(RP), public :: TIME_START              ! Start time [sec]
  real(RP), public :: TIME_END                ! End   time [sec]
  integer,  public :: TIME_NSTART             ! Time step at the start
  integer,  public :: TIME_NEND               ! Time step at the end

  real(RP), public :: TIME_CTIME              ! Current time [sec]
  integer,  public :: TIME_CSTEP              ! Current time step

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private :: TIME_STARTDATE(6) = (/ 0, 1, 1, 0, 0, 0 /)
  real(DP), private :: TIME_STARTMS      = 0.0_DP !< [millisec]
  integer,  private :: TIME_STARTDAY
  real(DP), private :: TIME_STARTSEC

  integer,  private :: TIME_ENDDATE(6)
  real(DP), private :: TIME_ENDMS
  integer,  private :: TIME_ENDDAY
  real(DP), private :: TIME_ENDSEC

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup the temporal scheme and time management
  subroutine TIME_setup
    use mod_adm, only: &
       ADM_CTL_FID, &
       ADM_proc_stop
    use scale_calendar, only: &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec, &
       CALENDAR_date2char
    use scale_time, only: &
       TIME_NOWDATE,     &
       TIME_NOWMS,       &
       TIME_NOWDAY,      &
       TIME_NOWSEC,      &
       TIME_NOWSTEP,     &
       TIME_NSTEP,       &
       TIME_OFFSET_YEAR
    implicit none

    character(len=ADM_NSYS) :: integ_type !--- integration method
    logical                 :: split      !--- time spliting flag
    real(RP)                 :: dtl        !--- delta t in large step
    integer                 :: lstep_max  !--- maximum number of large steps
    integer                 :: sstep_max  !--- division number in large step

    integer :: start_date(6) !< start date
    integer :: start_year    !< start year
    integer :: start_month   !< start month
    integer :: start_day     !< start day
    integer :: start_hour    !< start hour
    integer :: start_min     !< start min
    integer :: start_sec     !< start sec

    namelist / TIMEPARAM / &
         integ_type,  &
         split,       &
         dtl,         &
         lstep_max,   &
         sstep_max,   &
         start_date,  &
         start_year,  &
         start_month, &
         start_day,   &
         start_hour,  &
         start_min,   &
         start_sec

    character(len=27) :: startchardate
    character(len=27) :: endchardate

    integer :: ierr
    !---------------------------------------------------------------------------

    integ_type = TIME_integ_type
    split      = TIME_split
    dtl        = TIME_dtl
    lstep_max  = TIME_lstep_max
    sstep_max  = -999

    start_date(:) = -999
    start_year    = 0
    start_month   = 1
    start_day     = 1
    start_hour    = 0
    start_min     = 0
    start_sec     = 0

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[time]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=TIMEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** TIMEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist TIMEPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist TIMEPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=TIMEPARAM)

    !--- rewrite
    TIME_integ_type = integ_type
    TIME_split      = split
    TIME_dtl        = dtl
    TIME_lstep_max  = lstep_max

    if ( sstep_max == -999 )  then
       write(ADM_LOG_FID,*) 'TIME_integ_type is ', trim(TIME_integ_type)
       select case(TIME_integ_type)
       case('RK2')
          TIME_sstep_max = 4
       case('RK3')
          TIME_sstep_max = 6
        case('RK4')
          TIME_sstep_max = 8
        case('TRCADV')
          TIME_sstep_max = 0
       case default
          write(*,*) 'xxx Invalid TIME_INTEG_TYPE! STOP.'
       endselect
       write(ADM_LOG_FID,*) 'TIME_sstep_max is automatically set to: ', TIME_sstep_max
    else
       TIME_sstep_max = sstep_max
    endif
    TIME_dts = TIME_dtl / max(real(TIME_sstep_max,kind=RP),1.0_RP)

    if ( start_date(1) == -999 ) start_date(1) = start_year
    if ( start_date(2) == -999 ) start_date(2) = start_month
    if ( start_date(3) == -999 ) start_date(3) = start_day
    if ( start_date(4) == -999 ) start_date(4) = start_hour
    if ( start_date(5) == -999 ) start_date(5) = start_min
    if ( start_date(6) == -999 ) start_date(6) = start_sec

    TIME_OFFSET_YEAR = 0
    TIME_NOWSTEP     = 0
    TIME_NSTEP       = TIME_DTL

    TIME_STARTDATE(:) = start_date(:)
    TIME_STARTMS      = 0

    call CALENDAR_date2daysec( TIME_STARTDAY,     & ! [OUT]
                               TIME_STARTSEC,     & ! [OUT]
                               TIME_STARTDATE(:), & ! [IN]
                               TIME_STARTMS,      & ! [IN]
                               TIME_OFFSET_YEAR   ) ! [IN]

    TIME_start = CALENDAR_combine_daysec( TIME_STARTDAY, TIME_STARTSEC )

    call CALENDAR_date2char( startchardate,     & ! [OUT]
                             TIME_STARTDATE(:), & ! [IN]
                             TIME_STARTMS,      & ! [IN]
                             TIME_OFFSET_YEAR   ) ! [IN]

    TIME_END    = TIME_START  + TIME_LSTEP_MAX * TIME_DTL
    TIME_NSTART = 0
    TIME_NEND   = TIME_NSTART + TIME_LSTEP_MAX

    TIME_CTIME  = TIME_START
    TIME_CSTEP  = TIME_NSTART

    TIME_NOWDATE(:)   = TIME_STARTDATE(:)
    TIME_NOWMS        = TIME_STARTMS
    TIME_NOWDAY       = TIME_STARTDAY
    TIME_NOWSEC       = TIME_STARTSEC

    TIME_ENDDAY = 0
    TIME_ENDSEC = TIME_END

    call CALENDAR_adjust_daysec( TIME_ENDDAY, TIME_ENDSEC ) ! [INOUT]

    call CALENDAR_daysec2date( TIME_ENDDATE(:), & ! [OUT]
                               TIME_ENDMS,      & ! [OUT]
                               TIME_ENDDAY,     & ! [IN]
                               TIME_ENDSEC,     & ! [IN]
                               TIME_OFFSET_YEAR ) ! [IN]

    call CALENDAR_date2char( endchardate,     & ! [OUT]
                             TIME_ENDDATE(:), & ! [IN]
                             TIME_ENDMS,      & ! [IN]
                             TIME_OFFSET_YEAR ) ! [IN]

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '====== Time management ======'
    write(ADM_LOG_FID,*) '--- Time integration scheme (large step): ', trim(TIME_integ_type)
    write(ADM_LOG_FID,*) '--- Time interval for large step        : ', TIME_DTL
    write(ADM_LOG_FID,*) '--- Time interval for small step        : ', TIME_DTS
    write(ADM_LOG_FID,*) '--- Max steps of large step             : ', TIME_LSTEP_MAX
    write(ADM_LOG_FID,*) '--- Max steps of small step             : ', TIME_SSTEP_MAX
    write(ADM_LOG_FID,*) '--- Start time (sec)                    : ', TIME_START
    write(ADM_LOG_FID,*) '--- End time   (sec)                    : ', TIME_END
    write(ADM_LOG_FID,*) '--- Start time (date)                   : ', startchardate
    write(ADM_LOG_FID,*) '--- End time   (date)                   : ', endchardate
    write(ADM_LOG_FID,*) '--- total integration time              : ', TIME_END - TIME_START
    write(ADM_LOG_FID,*) '--- Time step at the start              : ', TIME_NSTART
    write(ADM_LOG_FID,*) '--- Time step at the end                : ', TIME_NEND

    return
  end subroutine TIME_setup

  !-----------------------------------------------------------------------------
  subroutine TIME_report
    use mod_adm, only: &
       ADM_prc_run_master, &
       ADM_prc_me
    use scale_calendar, only: &
       CALENDAR_date2char
    use scale_time, only: &
       TIME_NOWDATE,     &
       TIME_NOWMS,       &
       TIME_NOWSTEP,     &
       TIME_NSTEP,       &
       TIME_OFFSET_YEAR
    implicit none

    character(len=27) :: nowchardate
    !---------------------------------------------------------------------------

    call CALENDAR_date2char( nowchardate,     & ! [OUT]
                             TIME_NOWDATE(:), & ! [IN]
                             TIME_NOWMS,      & ! [IN]
                             TIME_OFFSET_YEAR ) ! [IN]

    write(ADM_LOG_FID,'(1x,3A,I6,A,I6)') '*** TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP
    if( ADM_prc_me == ADM_prc_run_master ) then
       write(*,'(1x,3A,I6,A,I6)') '*** TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP
    endif

    return
  end subroutine TIME_report

  !-----------------------------------------------------------------------------
  subroutine TIME_advance
    use scale_calendar, only: &
       CALENDAR_daysec2date,    &
       CALENDAR_adjust_daysec
    use scale_time, only: &
       TIME_NOWDATE,     &
       TIME_NOWMS,       &
       TIME_NOWDAY,      &
       TIME_NOWSEC,      &
       TIME_NOWSTEP,     &
       TIME_NSTEP,       &
       TIME_OFFSET_YEAR
    implicit none
    !---------------------------------------------------------------------------

    ! time advance
    TIME_CTIME = TIME_CTIME + TIME_DTL
    TIME_CSTEP = TIME_CSTEP + 1

    TIME_NOWSTEP = TIME_CSTEP
    TIME_NOWSEC  = TIME_CTIME

    ! reallocate day & sub-day
    call CALENDAR_adjust_daysec( TIME_NOWDAY, TIME_NOWSEC ) ! [INOUT]

    call CALENDAR_daysec2date( TIME_NOWDATE(:), & ! [OUT]
                               TIME_NOWMS,      & ! [OUT]
                               TIME_NOWDAY,     & ! [IN]
                               TIME_NOWSEC,     & ! [IN]
                               TIME_OFFSET_YEAR ) ! [IN]

    call TIME_report

    return
  end subroutine TIME_advance

end module mod_time
