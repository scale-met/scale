!-------------------------------------------------------------------------------
!> module ADMIN TIME
!!
!! @par Description
!!          Time management for SCALE-GM
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_time
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
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
  public :: TIME_setup
  public :: TIME_report
  public :: TIME_advance

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: TIME_integ_type = 'UNDEF'  ! Integration method in large steps
                                                    ! 'RK2'    ! Runge-Kutta 2nd
                                                    ! 'RK3'    ! Runge-Kutta 3rd
                                                    ! 'RK4'    ! Runge-Kutta 4th
                                                    ! 'TRCADV' ! Tracer advection only

  logical,  public :: TIME_split     = .true. ! Horizontally splitting?

  integer,  public :: TIME_lstep_max = 10     ! Max steps of large step
  integer,  public :: TIME_sstep_max          ! Max steps of small step

  real(DP), public :: TIME_dtl       = 5.0_RP ! Time interval for large step [sec]
  real(DP), public :: TIME_dts                ! Time interval for small step [sec]
  !
  real(DP), public :: TIME_START              ! Start time [sec]
  real(DP), public :: TIME_END                ! End   time [sec]

  real(DP), public :: TIME_CTIME              ! Current time [sec]
  integer,  public :: TIME_CSTEP              ! Current time step

  logical,  public :: TIME_DOend                !< finish program in this step?

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

  real(DP), private :: TIME_WALLCLOCK_START             ! Start time of wall clock             [sec]
  real(DP), private :: TIME_WALLCLOCK_LIMIT   = -1.0_DP ! Elapse time limit of wall clock time [sec]
  real(DP), private :: TIME_WALLCLOCK_SAFE    =  0.9_DP ! Safety coefficient for elapse time limit

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup the temporal scheme and time management
  subroutine TIME_setup
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_MPItime
    use scale_const, only: &
       UNDEF8 => CONST_UNDEF8
    use scale_calendar, only: &
       CALENDAR_date2daysec,    &
       CALENDAR_daysec2date,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec, &
       CALENDAR_unit2sec,       &
       CALENDAR_CFunits2sec,    &
       CALENDAR_date2char
    use scale_time, only: &
       TIME_DTSEC,                 &
       TIME_NOWDATE,               &
       TIME_NOWMS,                 &
       TIME_NOWDAY,                &
       TIME_NOWSEC,                &
       TIME_NOWDAYSEC,             &
       TIME_NOWSTEP,               &
       TIME_NSTEP,                 &
       TIME_DTSEC_WALLCLOCK_CHECK, &
       TIME_DSTEP_WALLCLOCK_CHECK, &
       TIME_OFFSET_YEAR,           &
       TIME_STARTDAYSEC
    implicit none

    character(len=H_SHORT) :: integ_type !--- integration method
    logical                :: split      !--- time spliting flag
    real(RP)               :: dtl        !--- delta t in large step
    integer                :: lstep_max  !--- maximum number of large steps
    integer                :: sstep_max  !--- division number in large step

    integer :: start_date(6) !< start date
    integer :: start_year    !< start year
    integer :: start_month   !< start month
    integer :: start_day     !< start day
    integer :: start_hour    !< start hour
    integer :: start_min     !< start min
    integer :: start_sec     !< start sec

    real(DP)               :: TIME_DT_WALLCLOCK_CHECK      = UNDEF8
    character(len=H_SHORT) :: TIME_DT_WALLCLOCK_CHECK_UNIT = ""

    NAMELIST / TIMEPARAM / &
       integ_type,                   &
       split,                        &
       dtl,                          &
       lstep_max,                    &
       sstep_max,                    &
       start_date,                   &
       start_year,                   &
       start_month,                  &
       start_day,                    &
       start_hour,                   &
       start_min,                    &
       start_sec,                    &
       TIME_DT_WALLCLOCK_CHECK,      &
       TIME_DT_WALLCLOCK_CHECK_UNIT, &
       TIME_WALLCLOCK_LIMIT,         &
       TIME_WALLCLOCK_SAFE

    real(DP)          :: TIME_DURATIONSEC
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
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[time]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=TIMEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** TIMEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist TIMEPARAM. STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist TIMEPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=TIMEPARAM)

    !--- rewrite
    TIME_integ_type = integ_type
    TIME_split      = split
    TIME_dtl        = dtl
    TIME_lstep_max  = lstep_max

    if ( sstep_max == -999 )  then
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
    else
       TIME_sstep_max = sstep_max
    endif
    TIME_dts = TIME_dtl / max(real(TIME_sstep_max,kind=DP),1.0_DP)

    if ( start_date(1) == -999 ) start_date(1) = start_year
    if ( start_date(2) == -999 ) start_date(2) = start_month
    if ( start_date(3) == -999 ) start_date(3) = start_day
    if ( start_date(4) == -999 ) start_date(4) = start_hour
    if ( start_date(5) == -999 ) start_date(5) = start_min
    if ( start_date(6) == -999 ) start_date(6) = start_sec

    TIME_STARTDATE(:) = start_date(:)
    TIME_STARTMS      = 0



    TIME_OFFSET_YEAR = TIME_STARTDATE(1)

    call CALENDAR_date2daysec( TIME_STARTDAY,     & ! [OUT]
                               TIME_STARTSEC,     & ! [OUT]
                               TIME_STARTDATE(:), & ! [IN]
                               TIME_STARTMS,      & ! [IN]
                               TIME_OFFSET_YEAR   ) ! [IN]

    call CALENDAR_date2char( startchardate,     & ! [OUT]
                             TIME_STARTDATE(:), & ! [IN]
                             TIME_STARTMS       ) ! [IN]

    TIME_STARTDAYSEC  = CALENDAR_combine_daysec( TIME_STARTDAY, TIME_STARTSEC )

    TIME_NOWDATE(:)   = TIME_STARTDATE(:)
    TIME_NOWMS        = TIME_STARTMS
    TIME_NOWDAY       = TIME_STARTDAY
    TIME_NOWSEC       = TIME_STARTSEC
    TIME_NOWDAYSEC    = CALENDAR_combine_daysec( TIME_NOWDAY, TIME_NOWSEC )

    TIME_ENDDAY       = TIME_STARTDAY

    TIME_DURATIONSEC  = TIME_lstep_max * TIME_dtl
    TIME_ENDSEC = TIME_STARTSEC + TIME_DURATIONSEC

    call CALENDAR_adjust_daysec( TIME_ENDDAY, TIME_ENDSEC ) ! [INOUT]

    call CALENDAR_daysec2date( TIME_ENDDATE(:), & ! [OUT]
                               TIME_ENDMS,      & ! [OUT]
                               TIME_ENDDAY,     & ! [IN]
                               TIME_ENDSEC,     & ! [IN]
                               TIME_OFFSET_YEAR ) ! [IN]

    call CALENDAR_date2char( endchardate,     & ! [OUT]
                             TIME_ENDDATE(:), & ! [IN]
                             TIME_ENDMS       ) ! [IN]

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Date/time setting ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)') '*** START Date     : ', startchardate
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)') '*** END   Date     : ', endchardate

    TIME_DTSEC   = TIME_dtl
    TIME_NSTEP   = TIME_lstep_max
    TIME_NOWSTEP = 1

    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** delta t (sec.) :', TIME_DTSEC
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I10)'  ) '*** No. of steps   :', TIME_NSTEP

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*)                     '*** Time interval for each processes (sec.)'
    if( IO_L ) write(IO_FID_LOG,*)                     '*** Atmosphere'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)')        '*** Dynamics small step (time)   : ', TIME_dts
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I10,A,I8,A)')   '***                     (step)   : ', TIME_sstep_max
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)')            '*** Dynamics large step (scheme) : ', trim(TIME_integ_type)

    TIME_CTIME  = TIME_NOWDAYSEC
    TIME_CSTEP  = TIME_NOWSTEP

    ! WALLCLOCK TERMINATOR SETUP
    if ( TIME_WALLCLOCK_LIMIT > 0.0_DP ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Wall clock time limit of execution is specified.'

       if ( TIME_DT_WALLCLOCK_CHECK == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_WALLCLOCK_CHECK. largest time step interval is used.'
          TIME_DTSEC_WALLCLOCK_CHECK = TIME_DTSEC
       else
          if ( TIME_DT_WALLCLOCK_CHECK_UNIT == '' ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_WALLCLOCK_CHECK_UNIT. SEC is used.'
             TIME_DT_WALLCLOCK_CHECK_UNIT = "SEC"
          endif
          call CALENDAR_unit2sec( TIME_DTSEC_WALLCLOCK_CHECK, TIME_DT_WALLCLOCK_CHECK, TIME_DT_WALLCLOCK_CHECK_UNIT )
          TIME_DTSEC_WALLCLOCK_CHECK = max( TIME_DTSEC_WALLCLOCK_CHECK, TIME_DTSEC )
       endif

       TIME_DSTEP_WALLCLOCK_CHECK = int( TIME_DTSEC_WALLCLOCK_CHECK / TIME_DTSEC )

       TIME_WALLCLOCK_SAFE  = max( min( TIME_WALLCLOCK_SAFE, 1.0_DP ), 0.0_DP )
       TIME_WALLCLOCK_START = PRC_MPItime()

       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.1,A)')      '*** This job stops after ', &
                                                          TIME_WALLCLOCK_LIMIT * TIME_WALLCLOCK_SAFE, ' seconds.'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Step interval for check     : ', TIME_DTSEC_WALLCLOCK_CHECK, &
                                                          ' (step interval=', TIME_DSTEP_WALLCLOCK_CHECK, ')'
    endif

    return
  end subroutine TIME_setup

  !-----------------------------------------------------------------------------
  subroutine TIME_report
    use scale_process, only: &
       PRC_MPItime
    use scale_calendar, only: &
       CALENDAR_date2char
    use scale_time, only: &
       TIME_NOWDATE, &
       TIME_NOWMS,   &
       TIME_NOWSTEP, &
       TIME_NSTEP
    implicit none

    real(DP)          :: WALLCLOCK_elapse
    character(len=27) :: nowchardate
    !---------------------------------------------------------------------------

    call CALENDAR_date2char( nowchardate,     & ! [OUT]
                             TIME_NOWDATE(:), & ! [IN]
                             TIME_NOWMS       ) ! [IN]

    if ( TIME_WALLCLOCK_LIMIT > 0.0_DP ) then
       WALLCLOCK_elapse = PRC_MPItime() - TIME_WALLCLOCK_START
       if( IO_L ) write(IO_FID_LOG,'(1x,3A,I7,A,I7,A,F10.1)') '*** TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP, &
                                                              ' WCLOCK:', WALLCLOCK_elapse
    else
       if( IO_L ) write(IO_FID_LOG,'(1x,3A,I7,A,I7)') '*** TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP
    endif

    return
  end subroutine TIME_report

  !-----------------------------------------------------------------------------
  subroutine TIME_advance( reverse )
    use scale_process, only: &
       PRC_IsMaster, &
       PRC_MPItime
    use scale_calendar, only: &
       CALENDAR_daysec2date,   &
       CALENDAR_adjust_daysec, &
       CALENDAR_combine_daysec
    use scale_time, only: &
       TIME_DTSEC,                &
       TIME_NOWDATE,              &
       TIME_NOWMS,                &
       TIME_NOWDAY,               &
       TIME_NOWSEC,               &
       TIME_NOWDAYSEC,            &
       TIME_NOWSTEP,              &
       TIME_NSTEP,                &
       TIME_OFFSET_YEAR,          &
       TIME_DSTEP_WALLCLOCK_CHECK
    implicit none

    logical, intent(in), optional :: reverse

    real(DP) :: WALLCLOCK_elapse
    integer  :: increment
    logical  :: exists
    !---------------------------------------------------------------------------

    TIME_DOend = .false.

    increment = 1
    if ( present(reverse) ) then
       if (reverse) then
          increment = - 1
       endif
    endif

    TIME_NOWSTEP = TIME_NOWSTEP + increment
    TIME_NOWDAY  = TIME_STARTDAY
    TIME_NOWSEC  = TIME_STARTSEC + real(TIME_NOWSTEP-1,kind=DP) * TIME_DTSEC

    ! reallocate day & sub-day
    call CALENDAR_adjust_daysec( TIME_NOWDAY, TIME_NOWSEC ) ! [INOUT]

    call CALENDAR_daysec2date( TIME_NOWDATE(:), & ! [OUT]
                               TIME_NOWMS,      & ! [OUT]
                               TIME_NOWDAY,     & ! [IN]
                               TIME_NOWSEC,     & ! [IN]
                               TIME_OFFSET_YEAR ) ! [IN]

    TIME_NOWDAYSEC = CALENDAR_combine_daysec( TIME_NOWDAY, TIME_NOWSEC )

    TIME_CTIME  = TIME_NOWDAYSEC
    TIME_CSTEP  = TIME_NOWSTEP

    if ( TIME_NOWSTEP > TIME_NSTEP ) then
       TIME_DOend = .true.
    endif

    if (       TIME_WALLCLOCK_LIMIT > 0.0_DP                     &
         .AND. mod(TIME_NOWSTEP-1,TIME_DSTEP_WALLCLOCK_CHECK) == 0 ) then
       WALLCLOCK_elapse = PRC_MPItime() - TIME_WALLCLOCK_START

       if ( WALLCLOCK_elapse > TIME_WALLCLOCK_LIMIT * TIME_WALLCLOCK_SAFE ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Elapse time limit is detected. Termination operation starts.'
          TIME_DOend = .true.
       endif
    endif

    ! QUIT file control
    if ( PRC_IsMaster ) then ! master node
       inquire(file='QUIT', exist=exists)

       if( exists ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** QUIT file is found. Termination operation starts.'
          TIME_DOend = .true.
       endif
    endif

    return
  end subroutine TIME_advance

end module mod_time
