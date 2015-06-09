!-------------------------------------------------------------------------------
!>
!! Time management module
!!
!! @par Description
!!         This module is for the time management.
!! @author  H.Tomita
!! @par History
!! @li      2004-02-17 (H.Tomita) Imported from igdc-4.33
!! @li      2004-05-31 (      )   Calculation of "num_of_iteration_[sl]step" are moved to mod[onestep].
!<
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

  logical, public :: TIME_SPLIT     = .true. ! Horizontally splitting?

  integer, public :: TIME_LSTEP_MAX = 10     ! Max steps of large step
  integer, public :: TIME_SSTEP_MAX          ! Max steps of small step

  real(RP), public :: TIME_DTL       = 5.0_RP   ! Time interval for large step [sec]
  real(RP), public :: TIME_DTS                ! Time interval for small step [sec]
  !
  real(RP), public :: TIME_START              ! Start time [sec]
  real(RP), public :: TIME_END                ! End   time [sec]
  integer, public :: TIME_NSTART             ! Time step at the start
  integer, public :: TIME_NEND               ! Time step at the end

  real(RP), public :: TIME_CTIME              ! Current time [sec]
  integer, public :: TIME_CSTEP              ! Current time step

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup the temporal scheme and time management
  subroutine TIME_setup
    use mod_adm, only: &
       ADM_CTL_FID, &
       ADM_proc_stop
    use mod_calendar, only: &
       calendar_yh2ss, &
       calendar_ss2cc
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

    character(len=20) :: HTIME_start
    character(len=20) :: HTIME_end

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
    call calendar_yh2ss( TIME_start, start_date )

    TIME_END    = TIME_START  + TIME_LSTEP_MAX * TIME_DTL
    TIME_NSTART = 0
    TIME_NEND   = TIME_NSTART + TIME_LSTEP_MAX

    TIME_CTIME  = TIME_START
    TIME_CSTEP  = TIME_NSTART

    !--- output the information for debug
    call calendar_ss2cc ( HTIME_start, TIME_START )
    call calendar_ss2cc ( HTIME_end,   TIME_END   )

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '====== Time management ======'
    write(ADM_LOG_FID,*) '--- Time integration scheme (large step): ', trim(TIME_integ_type)
    write(ADM_LOG_FID,*) '--- Time interval for large step        : ', TIME_DTL
    write(ADM_LOG_FID,*) '--- Time interval for small step        : ', TIME_DTS
    write(ADM_LOG_FID,*) '--- Max steps of large step             : ', TIME_LSTEP_MAX
    write(ADM_LOG_FID,*) '--- Max steps of small step             : ', TIME_SSTEP_MAX
    write(ADM_LOG_FID,*) '--- Start time (sec)                    : ', TIME_START
    write(ADM_LOG_FID,*) '--- End time   (sec)                    : ', TIME_END
    write(ADM_LOG_FID,*) '--- Start time (date)                   : ', HTIME_start
    write(ADM_LOG_FID,*) '--- End time   (date)                   : ', HTIME_end
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
    use mod_calendar, only: &
       calendar_ss2cc
    implicit none

    character(len=20) :: HTIME
    !---------------------------------------------------------------------------

    call calendar_ss2cc ( HTIME, TIME_CTIME )

    write(ADM_LOG_FID,*) '### TIME =', HTIME,'( step = ', TIME_CSTEP, ' )'
    if( ADM_prc_me == ADM_prc_run_master ) then
       write(*,*) '### TIME = ', HTIME,'( step = ', TIME_CSTEP, ' )'
    endif

    return
  end subroutine TIME_report

  !-----------------------------------------------------------------------------
  subroutine TIME_advance
    implicit none
    !---------------------------------------------------------------------------

    ! time advance
    TIME_CTIME = TIME_CTIME + TIME_DTL
    TIME_CSTEP = TIME_CSTEP + 1

    call TIME_report

    return
  end subroutine TIME_advance

end module mod_time
!-------------------------------------------------------------------------------
