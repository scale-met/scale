!-------------------------------------------------------------------------------
!> Module SNO (RM) PLUGIN
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          SCALE NetCDF Operator (SNO)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_snoplugin_timeave
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use mod_sno_h, only: &
     iteminfo
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SNOPLGIN_timeave_setup
  public :: SNOPLGIN_timeave_alloc
  public :: SNOPLGIN_timeave_dealloc
  public :: SNOPLGIN_timeave_store

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: SNOPLGIN_timeave_type        = 'OFF'  ! type of average
                                                                  ! 'OFF'    : disable
                                                                  ! 'NUMBER' : average by number
                                                                  ! 'DAILY'
                                                                  ! 'MONTHLY'
                                                                  ! 'ANNUAL'
  integer,                private :: SNOPLGIN_timeave_interval    = -1      ! number of averaging steps
  logical,                private :: SNOPLGIN_timeave_outorigdata = .false. ! output original (non-averaged) data ?

  integer,  private              :: timeave_counter    ! summation counter
  real(DP), private              :: timeave_time_start
  real(DP), private              :: timeave_time_end

  real(DP), private              :: timeave_startsec   ! from time_units
  integer,  private              :: timeave_refdate(6) ! reference date for Annual/Monthly/Daily average

  real(RP), private, allocatable :: counter_1d(:)
  real(RP), private, allocatable :: counter_2d(:,:)
  real(RP), private, allocatable :: counter_3d(:,:,:)

  type(iteminfo) :: avgdinfo

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_timeave_setup( &
       enable_plugin, &
       do_output      )
    use scale_prc, only: &
       PRC_abort
    implicit none

    logical, intent(out)   :: enable_plugin
    logical, intent(inout) :: do_output

    namelist / PARAM_SNOPLGIN_TIMEAVE / &
       SNOPLGIN_timeave_type,       &
       SNOPLGIN_timeave_interval,   &
       SNOPLGIN_timeave_outorigdata

    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SNOPLGIN_timeave_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SNOPLGIN_TIMEAVE,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
       LOG_INFO("SNOPLGIN_timeave_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SNOPLGIN_timeave_setup",*) 'Not appropriate names in namelist PARAM_SNOPLGIN_TIMEAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_SNOPLGIN_TIMEAVE)

    LOG_NEWLINE
    select case(SNOPLGIN_timeave_type)
    case('OFF')

       LOG_INFO("SNOPLGIN_timeave_setup",*) 'SNOPLGIN_timeave_type     : OFF'
       enable_plugin = .false.

    case('NUMBER')

       LOG_INFO("SNOPLGIN_timeave_setup",*) 'SNOPLGIN_timeave_type     : by number'
       enable_plugin = .true.

       if ( SNOPLGIN_timeave_interval < 1 ) then
          LOG_ERROR("SNOPLGIN_timeave_setup",*) 'Please set positive number for SNOPLGIN_timeave_interval : ', SNOPLGIN_timeave_interval
          call PRC_abort
       else
          LOG_INFO("SNOPLGIN_timeave_setup",*) 'SNOPLGIN_timeave_interval : ', SNOPLGIN_timeave_interval
       endif

    case('DAILY')

       LOG_INFO("SNOPLGIN_timeave_setup",*) 'SNOPLGIN_timeave_type     : Daily mean'
       enable_plugin = .true.

    case('MONTHLY')

       LOG_INFO("SNOPLGIN_timeave_setup",*) 'SNOPLGIN_timeave_type     : Monthly mean'
       enable_plugin = .true.

    case('ANNUAL')

       LOG_INFO("SNOPLGIN_timeave_setup",*) 'SNOPLGIN_timeave_type     : Annual mean'
       enable_plugin = .true.

    case default
       LOG_ERROR("SNOPLGIN_timeave_setup",*) 'the name of SNOPLGIN_timeave_type is not appropriate : ', trim(SNOPLGIN_timeave_type)
       LOG_ERROR_CONT(*) 'you can choose OFF,NUMBER,DAILY,MONTHLY,ANNUAL'
       call PRC_abort
    end select

    if ( enable_plugin ) then
       LOG_INFO("SNOPLGIN_timeave_setup",*) 'output original (non-averaged) data? : ', SNOPLGIN_timeave_outorigdata
       do_output = SNOPLGIN_timeave_outorigdata
    endif

    return
  end subroutine SNOPLGIN_timeave_setup

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_timeave_alloc( &
       dinfo, &
       debug  )
    use scale_calendar, only: &
       CALENDAR_daysec2date,   &
       CALENDAR_adjust_daysec, &
       CALENDAR_CFunits2sec
    use mod_sno_h, only: &
       iteminfo
    implicit none

    type(iteminfo), intent(in)  :: dinfo ! variable information               (input)
    logical,        intent(in)  :: debug

    integer  :: gout1, gout2, gout3

    integer  :: now_day     !< absolute day    (time=t)
    real(DP) :: now_sec     !< absolute second (time=t)
    integer  :: now_date(6) !< date            (time=t)
    real(DP) :: now_ms      !< subsecond       (time=t)

    integer  :: d
    !---------------------------------------------------------------------------

    if ( debug ) then
       LOG_INFO("SNOPLGIN_timeave_alloc",*) '[SNOPLGIN_timeave_alloc] allocate temporal array'
    endif

    if ( dinfo%dim_rank == 1 ) then

       gout1 = size(dinfo%VAR_1d(:))

       allocate( avgdinfo%VAR_1d(gout1) )
       allocate( counter_1d     (gout1) )
       avgdinfo%VAR_1d(:) = 0.0_RP
       counter_1d     (:) = 0.0_RP

    elseif( dinfo%dim_rank == 2 ) then

       gout1 = size(dinfo%VAR_2d(:,:),1)
       gout2 = size(dinfo%VAR_2d(:,:),2)

       allocate( avgdinfo%VAR_2d(gout1,gout2) )
       allocate( counter_2d     (gout1,gout2) )
       avgdinfo%VAR_2d(:,:) = 0.0_RP
       counter_2d     (:,:) = 0.0_RP

    elseif( dinfo%dim_rank == 3 ) then

       gout1 = size(dinfo%VAR_3d(:,:,:),1)
       gout2 = size(dinfo%VAR_3d(:,:,:),2)
       gout3 = size(dinfo%VAR_3d(:,:,:),3)

       allocate( avgdinfo%VAR_3d(gout1,gout2,gout3) )
       allocate( counter_3d     (gout1,gout2,gout3) )
       avgdinfo%VAR_3d(:,:,:) = 0.0_RP
       counter_3d     (:,:,:) = 0.0_RP

    endif



    ! reset step counter and reference date
    timeave_counter = 0

    timeave_startsec = CALENDAR_CFunits2sec( cftime=0.0_DP, cfunits=dinfo%time_units, offset_year=0 )

    now_day = 0
    now_sec = timeave_startsec + dinfo%time_end(1)

    call CALENDAR_adjust_daysec( now_day, now_sec ) ! [INOUT]

    call CALENDAR_daysec2date  ( now_date, now_ms,                & ! [OUT]
                                 now_day,  now_sec, offset_year=0 ) ! [IN]

    timeave_refdate(1) = now_date(1)
    timeave_refdate(2) = now_date(2)
    timeave_refdate(3) = now_date(3)
    timeave_refdate(4) = 0
    timeave_refdate(5) = 0
    timeave_refdate(6) = 0



    ! set variable information
    avgdinfo%varname     = dinfo%varname
    avgdinfo%description = dinfo%description
    avgdinfo%units       = dinfo%units
    avgdinfo%datatype    = dinfo%datatype
    avgdinfo%dim_rank    = dinfo%dim_rank
    avgdinfo%dim_name(:) = ""
    avgdinfo%dim_size(:) = 0
    do d = 1, dinfo%dim_rank
       avgdinfo%dim_name(d) = dinfo%dim_name(d)
       avgdinfo%dim_size(d) = dinfo%dim_size(d)
    enddo
    avgdinfo%transpose   = dinfo%transpose
    avgdinfo%time_units  = dinfo%time_units

    ! reset time information
    avgdinfo%step_nmax     = 0
    avgdinfo%time_start(:) = 0.0_DP
    avgdinfo%time_end  (:) = 0.0_DP
    avgdinfo%dt            = 0.0_DP

    return
  end subroutine SNOPLGIN_timeave_alloc

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_timeave_dealloc( &
       debug  )
    use mod_sno_h, only: &
       iteminfo
    implicit none

    logical,        intent(in)  :: debug
    !---------------------------------------------------------------------------

    if ( debug ) then
       LOG_INFO("SNOPLGIN_timeave_dealloc",*) '[SNOPLGIN_timeave_dealloc] deallocate temporal array'
    endif

    if( allocated(avgdinfo%VAR_1d) ) deallocate( avgdinfo%VAR_1d )
    if( allocated(avgdinfo%VAR_2d) ) deallocate( avgdinfo%VAR_2d )
    if( allocated(avgdinfo%VAR_3d) ) deallocate( avgdinfo%VAR_3d )
    if( allocated(counter_1d     ) ) deallocate( counter_1d      )
    if( allocated(counter_2d     ) ) deallocate( counter_2d      )
    if( allocated(counter_3d     ) ) deallocate( counter_3d      )

    return
  end subroutine SNOPLGIN_timeave_dealloc

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_timeave_store( &
       dirpath,       &
       basename,      &
       output_grads,  &
       nowrank,       &
       nowstep,       &
       nprocs_x_out,  &
       nprocs_y_out,  &
       nhalos_x,      &
       nhalos_y,      &
       hinfo,         &
       naxis,         &
       ainfo,         &
       dinfo,         &
       debug          )
    use scale_const, only: &
       CONST_UNDEF, &
       CONST_EPS
    use scale_calendar, only: &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec, &
       CALENDAR_CFunits2sec
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo,   &
       iteminfo
    use mod_sno_vars, only: &
       SNO_vars_write
    implicit none

    character(len=*), intent(in)    :: dirpath                               ! directory path                     (output)
    character(len=*), intent(in)    :: basename                              ! basename of file                   (output)
    logical,          intent(in)    :: output_grads
    integer,          intent(in)    :: nowrank                               ! current rank                       (output)
    integer,          intent(in)    :: nowstep                               ! current step                       (output)
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (output)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (output)
    integer,          intent(in)    :: nhalos_x                              ! number of x-axis halo grids        (global domain)
    integer,          intent(in)    :: nhalos_y                              ! number of y-axis halo grids        (global domain)
    type(commoninfo), intent(in)    :: hinfo                                 ! common information                 (input)
    integer,          intent(in)    :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(in)    :: ainfo(naxis)                          ! axis information                   (input)
    type(iteminfo),   intent(in)    :: dinfo                                 ! variable information               (input)
    logical,          intent(in)    :: debug


    integer  :: gout1, gout2, gout3

    integer  :: now_day      !< absolute day    (time=t)
    real(DP) :: now_sec      !< absolute second (time=t)
    integer  :: now_date(6)  !< date            (time=t)
    real(DP) :: now_ms       !< subsecond       (time=t)
    integer  :: next_day     !< absolute day    (time=t+1)
    real(DP) :: next_sec     !< absolute second (time=t+1)
    integer  :: next_date(6) !< date            (time=t+1)
    real(DP) :: next_ms      !< subsecond       (time=t+1)

    logical  :: do_output, finalize, add_rm_attr
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    ! summation

    if ( dinfo%dim_rank == 1 ) then

       gout1 = size(dinfo%VAR_1d(:))

       do i = 1, gout1
          if ( abs( dinfo%VAR_1d(i)-CONST_UNDEF ) > CONST_EPS ) then
             avgdinfo%VAR_1d(i) = avgdinfo%VAR_1d(i) + dinfo%VAR_1d(i)
             counter_1d     (i) = counter_1d     (i) + 1.0_RP
          endif
       enddo

    elseif( dinfo%dim_rank == 2 ) then

       gout1 = size(dinfo%VAR_2d(:,:),1)
       gout2 = size(dinfo%VAR_2d(:,:),2)

       do j = 1, gout2
       do i = 1, gout1
          if ( abs( dinfo%VAR_2d(i,j)-CONST_UNDEF ) > CONST_EPS ) then
             avgdinfo%VAR_2d(i,j) = avgdinfo%VAR_2d(i,j) + dinfo%VAR_2d(i,j)
             counter_2d     (i,j) = counter_2d     (i,j) + 1.0_RP
          endif
       enddo
       enddo

    elseif( dinfo%dim_rank == 3 ) then

       gout1 = size(dinfo%VAR_3d(:,:,:),1)
       gout2 = size(dinfo%VAR_3d(:,:,:),2)
       gout3 = size(dinfo%VAR_3d(:,:,:),3)

       do j = 1, gout3
       do i = 1, gout2
       do k = 1, gout1
          if ( abs( dinfo%VAR_3d(k,i,j)-CONST_UNDEF ) > CONST_EPS ) then
             avgdinfo%VAR_3d(k,i,j) = avgdinfo%VAR_3d(k,i,j) + dinfo%VAR_3d(k,i,j)
             counter_3d     (k,i,j) = counter_3d     (k,i,j) + 1.0_RP
          endif
       enddo
       enddo
       enddo

    endif

    ! check time to purge

    do_output = .false.

    if ( nowstep < dinfo%step_nmax ) then
       next_day = 0
       next_sec = CALENDAR_CFunits2sec( cftime=0.0_DP, cfunits=dinfo%time_units, offset_year=0 ) &
                + dinfo%time_end(nowstep+1)

       call CALENDAR_adjust_daysec( next_day, next_sec ) ! [INOUT]

       call CALENDAR_daysec2date  ( next_date, next_ms,                & ! [OUT]
                                    next_day,  next_sec, offset_year=0 ) ! [IN]

       next_date(4) = 0
       next_date(5) = 0
       next_date(6) = 0
    else
       next_date(:) = -1
    endif

    select case(SNOPLGIN_timeave_type)
    case('NUMBER')

       timeave_counter = timeave_counter + 1

       if ( timeave_counter == 1 ) then
          timeave_time_start = dinfo%time_start(nowstep)
       endif

       if ( timeave_counter == SNOPLGIN_timeave_interval ) then
          do_output = .true.
       endif

       if ( do_output ) then
          timeave_time_end   = dinfo%time_end  (nowstep)

          if ( dinfo%time_start(nowstep) == dinfo%time_end(nowstep) ) then ! snapshot
             avgdinfo%dt = ( timeave_time_end - timeave_time_start ) &
                         * real(SNOPLGIN_timeave_interval  ,kind=DP) &
                         / real(SNOPLGIN_timeave_interval-1,kind=DP)

             timeave_time_start = 0.5_DP * ( timeave_time_start + timeave_time_end )
             timeave_time_end   = timeave_time_start
          else ! average
             avgdinfo%dt = timeave_time_end - timeave_time_start
          endif
       endif

    case('DAILY')

       timeave_counter = timeave_counter + 1

       if ( timeave_counter == 1 ) then
          now_ms      = 0.0_RP
          now_date(:) = timeave_refdate(:)
          now_date(3) = now_date(3) - 1

          call CALENDAR_date2daysec( now_day,  now_sec,               & ! [OUT]
                                     now_date, now_ms,  offset_year=0 ) ! [IN]

          timeave_time_start = CALENDAR_combine_daysec( now_day, now_sec ) - timeave_startsec
       endif

       if (      next_date(1) /= timeave_refdate(1) &
            .OR. next_date(2) /= timeave_refdate(2) &
            .OR. next_date(3) /= timeave_refdate(3) ) then
          do_output = .true.
       endif

       if ( do_output ) then
          now_ms      = 0.0_RP
          now_date(:) = timeave_refdate(:)

          call CALENDAR_date2daysec( now_day,  now_sec,               & ! [OUT]
                                     now_date, now_ms,  offset_year=0 ) ! [IN]

          timeave_time_end = CALENDAR_combine_daysec( now_day, now_sec ) - timeave_startsec

          avgdinfo%dt = timeave_time_end - timeave_time_start
       endif

    case('MONTHLY')

       timeave_counter = timeave_counter + 1

       if ( timeave_counter == 1 ) then
          now_ms      = 0.0_RP
          now_date(:) = timeave_refdate(:)
          now_date(2) = now_date(2) - 1
          now_date(3) = 1

          call CALENDAR_date2daysec( now_day,  now_sec,               & ! [OUT]
                                     now_date, now_ms,  offset_year=0 ) ! [IN]

          timeave_time_start = CALENDAR_combine_daysec( now_day, now_sec ) - timeave_startsec
       endif

       if (      next_date(1) /= timeave_refdate(1) &
            .OR. next_date(2) /= timeave_refdate(2) ) then
          do_output = .true.
       endif

       if ( do_output ) then
          now_ms      = 0.0_RP
          now_date(:) = timeave_refdate(:)
          now_date(3) = 1

          call CALENDAR_date2daysec( now_day,  now_sec,               & ! [OUT]
                                     now_date, now_ms,  offset_year=0 ) ! [IN]

          timeave_time_end = CALENDAR_combine_daysec( now_day, now_sec ) - timeave_startsec

          avgdinfo%dt = timeave_time_end - timeave_time_start
       endif

    case('ANNUAL')

       timeave_counter = timeave_counter + 1

       if ( timeave_counter == 1 ) then
          now_ms      = 0.0_RP
          now_date(:) = timeave_refdate(:)
          now_date(1) = now_date(1) - 1
          now_date(2) = 1
          now_date(3) = 1

          call CALENDAR_date2daysec( now_day,  now_sec,               & ! [OUT]
                                     now_date, now_ms,  offset_year=0 ) ! [IN]

          timeave_time_start = CALENDAR_combine_daysec( now_day, now_sec ) - timeave_startsec
       endif

       if (      next_date(1) /= timeave_refdate(1) ) then
          do_output = .true.
       endif

       if ( do_output ) then
          now_ms      = 0.0_RP
          now_date(:) = timeave_refdate(:)
          now_date(2) = 1
          now_date(3) = 1

          call CALENDAR_date2daysec( now_day,  now_sec,               & ! [OUT]
                                     now_date, now_ms,  offset_year=0 ) ! [IN]

          timeave_time_end = CALENDAR_combine_daysec( now_day, now_sec ) - timeave_startsec

          avgdinfo%dt = timeave_time_end - timeave_time_start
       endif

    end select

    ! output

    if ( do_output ) then
       ! store time information
       avgdinfo%step_nmax = avgdinfo%step_nmax + 1

       avgdinfo%time_start(avgdinfo%step_nmax) = timeave_time_start
       avgdinfo%time_end  (avgdinfo%step_nmax) = timeave_time_end

       ! store data
       if ( dinfo%dim_rank == 1 ) then

          gout1 = size(dinfo%VAR_1d(:))

          do i = 1, gout1
             if ( counter_1d(i) > 0.0_RP ) then
                avgdinfo%VAR_1d(i) = avgdinfo%VAR_1d(i) / counter_1d(i)
             endif
          enddo

       elseif( dinfo%dim_rank == 2 ) then

          gout1 = size(dinfo%VAR_2d(:,:),1)
          gout2 = size(dinfo%VAR_2d(:,:),2)

          do j = 1, gout2
          do i = 1, gout1
             if ( counter_2d(i,j) > 0.0_RP ) then
                avgdinfo%VAR_2d(i,j) = avgdinfo%VAR_2d(i,j) / counter_2d(i,j)
             endif
          enddo
          enddo

       elseif( dinfo%dim_rank == 3 ) then

          gout1 = size(dinfo%VAR_3d(:,:,:),1)
          gout2 = size(dinfo%VAR_3d(:,:,:),2)
          gout3 = size(dinfo%VAR_3d(:,:,:),3)

          do j = 1, gout3
          do i = 1, gout2
          do k = 1, gout1
             if ( counter_3d(k,i,j) > 0.0_RP ) then
                avgdinfo%VAR_3d(k,i,j) = avgdinfo%VAR_3d(k,i,j) / counter_3d(k,i,j)
             endif
          enddo
          enddo
          enddo

       endif

       LOG_INFO("SNOPLGIN_timeave_store",'(1x,A,I6)') '++ output tave = ', avgdinfo%step_nmax

       finalize    = ( nowstep == dinfo%step_nmax )
       add_rm_attr = .true.

       call SNO_vars_write( dirpath,                    & ! [IN] from namelist
                            basename,                   & ! [IN] from namelist
                            output_grads,               & ! [IN] from namelist
                            nowrank,                    & ! [IN]
                            avgdinfo%step_nmax,         & ! [IN]
                            finalize,                   & ! [IN]
                            add_rm_attr,                & ! [IN]
                            nprocs_x_out, nprocs_y_out, & ! [IN] from namelist
                            nhalos_x,     nhalos_y,     & ! [IN] from SNO_file_getinfo
                            hinfo,                      & ! [IN] from SNO_file_getinfo
                            naxis,                      & ! [IN] from SNO_file_getinfo
                            ainfo(:),                   & ! [IN] from SNO_axis_getinfo
                            avgdinfo,                   & ! [IN] from SNO_vars_getinfo
                            debug                       ) ! [IN]

       ! reset data buffer
       if ( dinfo%dim_rank == 1 ) then
          avgdinfo%VAR_1d(:) = 0.0_RP
          counter_1d     (:) = 0.0_RP
       elseif( dinfo%dim_rank == 2 ) then
          avgdinfo%VAR_2d(:,:) = 0.0_RP
          counter_2d     (:,:) = 0.0_RP
       elseif( dinfo%dim_rank == 3 ) then
          avgdinfo%VAR_3d(:,:,:) = 0.0_RP
          counter_3d     (:,:,:) = 0.0_RP
       endif

       ! reset step counter and reference date
       timeave_counter    = 0
       timeave_refdate(:) = next_date(:)
    endif

    return
  end subroutine SNOPLGIN_timeave_store

end module mod_snoplugin_timeave
