!-------------------------------------------------------------------------------
!> Module external data
!!
!! @par Description
!!          General module for reading external-data
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_extdata
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
  !++ Public procedures
  !
  public :: extdata_setup
  public :: extdata_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: data_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private            :: max_extdata            !--- max number of external data
  integer, private, parameter :: max_num_of_data = 2500 !--- max time step num

  !--- type definition of information of external data
  type, private :: extdatainfo
     character(len=H_LONG)  :: fname             !--- data file name
     character(len=H_LONG)  :: dataname          !--- data name
     character(len=H_SHORT) :: input_io_mode     !--- io mode
     integer                :: input_size        !--- double/single precision
     character(len=H_SHORT) :: layer_type        !--- type of layer : 'ATM' or 'SFC'
     character(len=H_SHORT) :: layername         !--- name of layer
     integer                :: kall              !--- number of layer
     integer                :: num_of_data       !--- number of data
     integer,  pointer      :: data_date(:,:)    !--- date of each data piece
     real(DP), pointer      :: data_time(:)      !--- time[sec] of each data piecce
     integer                :: data_rec(2)       !--- data record ( forward & backward )
     integer                :: fix_rec           !--- record number if it's fixed
     logical                :: opt_fix_rec       !--- flag for fix record number
     logical                :: opt_monthly_cnst  !--- no time interpolation in a month
     logical                :: opt_periodic_year !--- flag for periodic year
     real(RP)               :: defval            !--- default value
     real(RP), pointer      :: v(:,:,:,:)        !--- data stoarege for regular region
     real(RP), pointer      :: v_pl(:,:,:,:)     !--- data stoarege for poler region
  end type extdatainfo

  type(extdatainfo), allocatable, private :: info(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine extdata_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_calendar, only: &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec
    use mod_adm, only: &
       ADM_KNONE,    &
       ADM_kall,     &
       ADM_gall,     &
       ADM_lall,     &
       ADM_gall_pl,  &
       ADM_lall_pl
    use mod_fio, only: &
       FIO_seek
    use mod_time, only: &
       ctime => TIME_CTIME
    implicit none

    character(len=H_LONG)  :: fname
    character(len=H_LONG)  :: dataname
    character(len=H_SHORT) :: input_io_mode
    integer                :: input_size
    character(len=H_SHORT) :: layer_type
    character(len=H_SHORT) :: layername
    integer                :: nlayer
    integer                :: num_of_data
    integer                :: data_date(6,max_num_of_data)
    integer                :: fix_rec
    logical                :: opt_fix_rec
    logical                :: opt_monthly_cnst
    logical                :: opt_periodic_year
    real(RP)               :: defval

    integer :: ddata_date(6)
    logical :: opt_increment_date

    namelist /nm_extdata/   &
         fname,             &
         dataname,          &
         input_io_mode,     &
         input_size,        &
         layer_type,        &
         layername,         &
         nlayer,            &
         num_of_data,       &
         data_date,         &
         ddata_date,        &
         opt_increment_date,&
         fix_rec,           &
         opt_fix_rec,       &
         opt_monthly_cnst,  &
         opt_periodic_year, &
         defval

    integer  :: cdate(6)

    integer  :: extday, offset_year
    real(DP) :: extsec, extms

    integer  :: ierr
    integer  :: np, im
    !---------------------------------------------------------------------------

    extday      = 0
    extsec      = ctime
    offset_year = 0

    call CALENDAR_adjust_daysec( extday, extsec ) ! [INOUT]

    call CALENDAR_daysec2date( cdate(:),   & ! [OUT]
                               extms,      & ! [OUT]
                               extday,     & ! [IN]
                               extsec,     & ! [IN]
                               offset_year ) ! [IN]

    !--- search the number of data files.
    max_extdata = 0
    rewind(IO_FID_CONF)
    do
       read(IO_FID_CONF, nml=nm_extdata, iostat=ierr)

       if ( ierr>0 ) then
          if( IO_L ) write(IO_FID_LOG,*) 'Msg : Sub[extdata_setup]/Mod[mod_extdata]'
          if( IO_L ) write(IO_FID_LOG,*) ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
          call PRC_MPIstop
       endif

       if(ierr < 0) exit

       max_extdata = max_extdata + 1
    enddo

    !--- allocation of data information array
    allocate( info(max_extdata) )

    !--- read namelist, again.
    rewind(IO_FID_CONF)
    do np=1, max_extdata

       !--- intialization
       fname             = ''
       dataname          = ''
       input_io_mode     = 'ADVANCED'
       input_size        = 8
       layer_type        = ''
       layername         = 'ZSSFC1'
       num_of_data       = 1
       data_date(:,:)    = 1
       fix_rec           = 1
       opt_fix_rec       = .false.
       opt_monthly_cnst  = .false.
       opt_periodic_year = .false.
       defval            = 0.0_RP
       ddata_date(:)     = 0
       opt_increment_date= .false.
       !
       !--- read namelist
       read(IO_FID_CONF, nml=nm_extdata, iostat=ierr)

       if ( trim(layer_type) == 'ATM' ) then
          info(np)%kall = ADM_kall
       elseif( trim(layer_type) == 'SFC' ) then
          info(np)%kall = ADM_KNONE
       elseif( trim(layer_type) == 'NUM' ) then
          info(np)%kall = nlayer
       else
          write(*,*) 'xxx [extdata] invlalid type of layer_type.', trim(layer_type)
          call PRC_MPIstop
       endif

       ! IO_seek verifies;
       !  info(np)%layername
       !  info(np)%kall
       !
       ! IO_seek overwrites;
       !  info(np)%input_size
       !  info(np)%num_of_data
       !  info(np)%data_date
       !  info(np)%data_rec(1)
       if ( input_io_mode == 'ADVANCED' ) then

          call FIO_seek( info(np)%data_rec(1), & ! [OUT]
                         num_of_data,          & ! [INOUT]
                         data_date,            & ! [INOUT]
                         input_size,           & ! [INOUT]
                         fname,                & ! [IN]
                         dataname,             & ! [IN]
                         layername,            & ! [IN]
                         1,                    & ! [IN]
                         info(np)%kall,        & ! [IN]
                         ctime,                & ! [IN]
                         cdate,                & ! [IN]
                         opt_periodic_year     ) ! [IN]

       else
          write(*,*) 'xxx [extdata] Invalid input_io_mode!', trim(input_io_mode)
          call PRC_MPIstop
       endif

       allocate( info(np)%data_date(6,num_of_data) )
       allocate( info(np)%data_time(  num_of_data) )
       allocate( info(np)%v   (ADM_gall,   info(np)%kall,ADM_lall   ,2) )
       allocate( info(np)%v_pl(ADM_gall_pl,info(np)%kall,ADM_lall_pl,2) )

       !--- store information
       info(np)%fname             = trim(fname)
       info(np)%dataname          = trim(dataname)
       info(np)%input_io_mode     = trim(input_io_mode)
       info(np)%input_size        = input_size
       info(np)%layer_type        = trim(layer_type)
       info(np)%layername         = trim(layername)
       info(np)%num_of_data       = num_of_data
       info(np)%fix_rec           = fix_rec
       info(np)%opt_fix_rec       = opt_fix_rec
       info(np)%opt_monthly_cnst  = opt_monthly_cnst
       info(np)%opt_periodic_year = opt_periodic_year
       info(np)%defval            = defval

       if ( num_of_data == 1 ) then
          info(np)%opt_fix_rec = .true.
          info(np)%fix_rec     = 1
       endif

       if ( opt_monthly_cnst .AND. num_of_data /= 12 ) then
          write(*,*) 'Msg : Sub[extdata_setup]/Mod[mod_extdata]'
          write(*,*) '---- ERROR! : inconsistent data num in monthly fix version!'
          call PRC_MPIstop
       endif

       info(np)%data_date(:,1:num_of_data) = data_date(:,1:num_of_data)
       do im = 1, num_of_data
          extms       = 0
          offset_year = 0

          call CALENDAR_date2daysec( extday,                   & ! [OUT]
                                     extsec,                   & ! [OUT]
                                     info(np)%data_date(:,im), & ! [IN]
                                     extms,                    & ! [IN]
                                     offset_year               ) ! [IN]

          info(np)%data_time(im) = CALENDAR_combine_daysec( extday, extsec )
       enddo

    enddo

    !---- set initial data record
    do np = 1, max_extdata

       if(info(np)%opt_fix_rec) then

          info(np)%data_rec(1) = info(np)%fix_rec
          info(np)%data_rec(2) = info(np)%fix_rec

       elseif(info(np)%opt_monthly_cnst) then

          info(np)%data_rec(1) = cdate(2)
          info(np)%data_rec(2) = cdate(2)

       elseif(info(np)%opt_periodic_year) then

          if( info(np)%data_rec(1) == 1 ) then
             info(np)%data_rec(2) = info(np)%num_of_data
          else
             info(np)%data_rec(2) = info(np)%data_rec(1)-1
          endif

       else !--- default

          if ( info(np)%data_rec(1) == 1 ) then
             write(*,*) 'xxx [extdata] data time is not consistent with the simulation time! : ', &
                        trim(info(np)%dataname )
             call PRC_MPIstop
          else !--- default
             info(np)%data_rec(2) = info(np)%data_rec(1)-1
          endif

       endif

       !--- read initial data
       call data_read( np )

    enddo

    !--- output information
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Msg : Sub[extdata_setup]/Mod[mod_extdata]'
    if( IO_L ) write(IO_FID_LOG,*) '===================================================='
    if( IO_L ) write(IO_FID_LOG,*) '--- Number of maximum external data  : ',max_extdata
    do np = 1, max_extdata
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,A)') '--- variable [',np,'] : ', trim(info(np)%dataname)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '===================================================='

    do np = 1, max_extdata
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A)') '============ External file NO. : ',np,' ============='
       if( IO_L ) write(IO_FID_LOG,*) '--- fname             : ', trim(info(np)%fname)
       if( IO_L ) write(IO_FID_LOG,*) '--- dataname          : ', trim(info(np)%dataname)
       if( IO_L ) write(IO_FID_LOG,*) '--- input_io_mode     : ', trim(info(np)%input_io_mode)
       if( IO_L ) write(IO_FID_LOG,*) '--- input_size        : ', info(np)%input_size
       if( IO_L ) write(IO_FID_LOG,*) '--- layer_type        : ', info(np)%layer_type
       if( IO_L ) write(IO_FID_LOG,*) '--- layername         : ', info(np)%layername
       if( IO_L ) write(IO_FID_LOG,*) '--- num_of_data       : ', info(np)%num_of_data
       do im = 1, info(np)%num_of_data
          if( IO_L ) write(IO_FID_LOG,'(1x,A,6(I4,1x))') '--- data_date         : ', info(np)%data_date(:,im)
       enddo
       if( IO_L ) write(IO_FID_LOG,*) '--- opt_fix_rec       : ', info(np)%opt_fix_rec
       if( IO_L ) write(IO_FID_LOG,*) '--- opt_monthly_cnst  : ', info(np)%opt_monthly_cnst
       if( IO_L ) write(IO_FID_LOG,*) '--- opt_periodic_year : ', info(np)%opt_periodic_year
       if( IO_L ) write(IO_FID_LOG,*) '--- defval            : ', info(np)%defval
       if ( info(np)%input_io_mode == 'ADVANCED' ) then
         if( IO_L ) write(IO_FID_LOG,*) '--- first step        : ', info(np)%data_rec(1)
       endif
       if( IO_L ) write(IO_FID_LOG,*) '====================================================='
    enddo

    return
  end subroutine extdata_setup

  !-----------------------------------------------------------------------------
  subroutine extdata_update( &
       gdata,    & !--- INOUT : data
       DNAME,    & !--- IN    : data name
       l_region, & !--- IN    : index of region
       ctime,    & !--- IN    : record number of data on current time
       eflag     )
    use scale_process, only: &
       PRC_MPIstop
    use scale_calendar, only: &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec
    use mod_adm, only: &
       ADM_gall_in, &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    real(RP),         intent(inout) :: gdata(:,:) ! data is inout to retain initilized value.
    character(len=*), intent(in)    :: DNAME      ! data name
    integer,          intent(in)    :: l_region
    real(DP),         intent(in)    :: ctime      ! current time
    logical,          intent(out)   :: eflag
    !
    real(RP) :: dt !-- delta t between two timestep data
    real(RP) :: wt !-- t weight of two timestep data

    !--- data ID
    integer  :: np
    integer  :: data_date_prev(6), cdate(6)
    real(DP) :: data_time_prev

    integer  :: extday, offset_year
    real(DP) :: extsec, extms

    integer  :: kall, gall
    integer  :: im, i, j, k, n
    !---------------------------------------------------------------------------

    eflag = .false.
    !--- searching the data ID
    np = 0
    do n = 1, max_extdata
       if ( trim(DNAME) == trim(info(n)%dataname) ) then
          np = n
          eflag = .true.
          exit
       endif
    enddo
    if ( np == 0 ) return

    !--- error handling for number of array
    gall = size(gdata,1)
    kall = size(gdata,2)
    if (      kall /= info(np)%kall &
         .OR. gall /= ADM_gall_in   ) then
       write(*,*) 'xxx Array size of gdata is not consistent!', trim(DNAME)
       call PRC_MPIstop
    endif

    !--- update external data
    if (info(np)%opt_fix_rec) then
       !--- no time-advance
    else
       !--- data time advance
       if ( info(np)%opt_monthly_cnst ) then

          extday      = 0
          extsec      = ctime
          offset_year = 0

          call CALENDAR_adjust_daysec( extday, extsec ) ! [INOUT]

          call CALENDAR_daysec2date( cdate(:),   & ! [OUT]
                                     extms,      & ! [OUT]
                                     extday,     & ! [IN]
                                     extsec,     & ! [IN]
                                     offset_year ) ! [IN]

          if (cdate(2) /= info(np)%data_rec(1)) then
             info(np)%data_rec(1) = cdate(2)
             info(np)%data_rec(2) = info(np)%data_rec(1)

             call data_read( np )
          endif

       else !--- default

          if ( ctime > info(np)%data_time(info(np)%data_rec(1)) ) then
             !<-- current time pass the current data

             if( IO_L ) write(IO_FID_LOG,*) '*** Update external data :',trim(info(np)%dataname)

             !--- increment of data_rec
             info(np)%data_rec(2) = info(np)%data_rec(1)
             info(np)%data_rec(1) = info(np)%data_rec(1) + 1
             if (       ( info(np)%data_rec(1) > info(np)%num_of_data ) &
                  .AND. ( .not. info(np)%opt_periodic_year )            ) then

                write(*,*) 'xxx [extdata] Current time exceeded the time range of the input data.'
                call PRC_MPIstop

             elseif( ( info(np)%data_rec(1) > info(np)%num_of_data ) .and. &
                     ( info(np)%opt_periodic_year ) ) then
                info(np)%data_rec(1) = info(np)%data_rec(1) - info(np)%num_of_data !-- 1
                info(np)%data_rec(2) = info(np)%num_of_data

                !--- re-setting of data_date and data_time
                info(np)%data_date(1,1:info(np)%num_of_data) = info(np)%data_date(1,1:info(np)%num_of_data) + 1
                do im = 1, info(np)%num_of_data
                   extms       = 0
                   offset_year = 0

                   call CALENDAR_date2daysec( extday,                   & ! [OUT]
                                              extsec,                   & ! [OUT]
                                              info(np)%data_date(:,im), & ! [IN]
                                              extms,                    & ! [IN]
                                              offset_year               ) ! [IN]

                   info(np)%data_time(im) = CALENDAR_combine_daysec( extday, extsec )
                enddo

                !--- output of message.
                if( IO_L ) write(IO_FID_LOG,*) '*** data date is updated as follows.'
                do im = 1,info(np)%num_of_data
                   if( IO_L ) write(IO_FID_LOG,*) ' ----- data_date(',im,') : ', info(np)%data_date(:,im)
                enddo

             endif

             call data_read( np )
          endif
       endif
    endif

    !--- store the raw data to gdata.
    if ( (info(np)%opt_fix_rec).or.(info(np)%opt_monthly_cnst) ) then
       wt = 1.0_RP
    elseif(info(np)%data_rec(1) == 1) then !<--- this case is only periodic one.
       data_date_prev(:) = info(np)%data_date(:,info(np)%num_of_data)
       data_date_prev(1) = data_date_prev(1) - 1

       extms       = 0
       offset_year = 0

       call CALENDAR_date2daysec( extday,            & ! [OUT]
                                  extsec,            & ! [OUT]
                                  data_date_prev(:), & ! [IN]
                                  extms,             & ! [IN]
                                  offset_year        ) ! [IN]

       data_time_prev = CALENDAR_combine_daysec( extday, extsec )

       dt = info(np)%data_time(info(np)%data_rec(1)) &
          - data_time_prev
       wt = ( ctime - data_time_prev ) / dt
    else !--- default
       dt = info(np)%data_time(info(np)%data_rec(1)) &
          - info(np)%data_time(info(np)%data_rec(2))
       wt = ( ctime - info(np)%data_time(info(np)%data_rec(2)) ) / dt
    endif

    do k = 1, kall
       n = 1
       do j = ADM_gmin, ADM_gmax+1
       do i = ADM_gmin, ADM_gmax+1
          gdata(n,k) = info(np)%v(suf(i,j),k,l_region,1) * (       wt) &
                     + info(np)%v(suf(i,j),k,l_region,2) * (1.0_RP-wt)

          n = n + 1
       enddo
       enddo
    enddo

    return
  end subroutine extdata_update

  !-----------------------------------------------------------------------------
  subroutine data_read( np )
    use mod_comm, only: &
       COMM_var
    use mod_fio, only: &
       FIO_input
    implicit none

    integer, intent(in) :: np

    integer :: n
    !---------------------------------------------------------------------------

    info(np)%v   (:,:,:,:) = info(np)%defval
    info(np)%v_pl(:,:,:,:) = info(np)%defval

    do n = 1, 2 !--- forward & backward
       if ( info(np)%input_io_mode == 'ADVANCED' ) then

          call FIO_input( info(np)%v(:,:,:,n), & ! [OUT]
                          info(np)%fname,      & ! [IN]
                          info(np)%dataname,   & ! [IN]
                          info(np)%layername,  & ! [IN]
                          1,info(np)%kall,     & ! [IN]
                          info(np)%data_rec(n) ) ! [IN]

       endif
    enddo

    call COMM_var( info(np)%v, info(np)%v_pl, info(np)%kall, 2 )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '===== READ EXTERNAL DATA ============================'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,A)') '--- variable [',np,'] : ', trim(info(np)%dataname)
    if( IO_L ) write(IO_FID_LOG,*) '--- forward  data step(record) number : ',info(np)%data_rec(1)
    if( IO_L ) write(IO_FID_LOG,*) '--- backward data step(record) number : ',info(np)%data_rec(2)
    if( IO_L ) write(IO_FID_LOG,*) '====================================================='

    return
  end subroutine data_read

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_extdata
