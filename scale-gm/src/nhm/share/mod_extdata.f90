!-------------------------------------------------------------------------------
!
!+  external data module
!
!-------------------------------------------------------------------------------
module mod_extdata
  !-----------------------------------------------------------------------------
  !
  !++ Description: General module of reading external-data.
  !                This modules was remade, based on the old modules[mod_extdata],
  !                [mod_o3var], and [mod_ocean_mixedlayer_exdata.f90]
  !                for integration of input routine of external files.
  !                Contributer : M.Satoh, T.Mitsui, Y.Niwa, W. Yanase, H. Tomita
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.01      09-01-22   H.Tomita  : Add this module.
  !                11-09-03   H.Yashiro : New I/O
  !                12-02-01   T.Seiki   : add option "increment date"
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_MAXFNAME, &
     ADM_NSYS
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
     character(ADM_MAXFNAME) :: fname             !--- data file name
     character(ADM_MAXFNAME) :: dataname          !--- data name
     character(ADM_MAXFNAME) :: input_io_mode     !--- io mode                  ! [add] H.Yashiro 20110826
     integer                 :: input_size        !--- double/single precision
     character(ADM_NSYS)     :: layer_type        !--- type of layer : 'ATM' or 'SFC'
     character(ADM_NSYS)     :: layername         !--- name of layer            ! [add] H.Yashiro 20110826
     integer                 :: kall              !--- number of layer
     integer                 :: num_of_data       !--- number of data
     integer, pointer        :: data_date(:,:)    !--- date of each data piece
     real(RP), pointer        :: data_time(:)      !--- time[sec] of each data piecce
     integer                 :: data_rec(2)       !--- data record ( forward & backward )
     integer                 :: fix_rec           !--- record number if it's fixed
     logical                 :: opt_fix_rec       !--- flag for fix record number
     logical                 :: opt_monthly_cnst  !--- no time interpolation in a month
     logical                 :: opt_periodic_year !--- flag for periodic year
     real(RP)                 :: defval            !--- default value
     real(RP), pointer        :: v(:,:,:,:)        !--- data stoarege for regular region
     real(RP), pointer        :: v_pl(:,:,:,:)     !--- data stoarege for poler region
  end type extdatainfo

  type(extdatainfo), allocatable, private :: info(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine extdata_setup
    use mod_adm, only : &
       ADM_CTL_FID,  &
       ADM_LOG_FID,  &
       ADM_KNONE,    &
       ADM_kall,     &
       ADM_gall,     &
       ADM_lall,     &
       ADM_gall_pl,  &
       ADM_lall_pl,  &
       ADM_proc_stop
    use mod_fio, only: &
       FIO_seek
    use mod_calendar, only: &
       calendar_yh2ss, &
       calendar_ss2yh
    use mod_time, only: &
       ctime => TIME_CTIME
    implicit none

    character(len=ADM_MAXFNAME) :: fname
    character(len=ADM_MAXFNAME) :: dataname
    character(len=ADM_MAXFNAME) :: input_io_mode ! [add] H.Yashiro 20110826
    integer                     :: input_size
    character(LEN=ADM_NSYS)     :: layer_type
    character(LEN=ADM_NSYS)     :: layername ! [add] H.Yashiro 20110826
    integer                     :: nlayer    ! [add] H.Yashiro 20131030
    integer                     :: num_of_data
    integer                     :: data_date(6,max_num_of_data)
    integer                     :: fix_rec
    logical                     :: opt_fix_rec
    logical                     :: opt_monthly_cnst
    logical                     :: opt_periodic_year
    real(RP)                     :: defval
    ! [Add] 12/02/01 T.Seiki
    integer :: ddata_date(6)
    logical :: opt_increment_date

    namelist /nm_extdata/   &
         fname,             &
         dataname,          &
         input_io_mode,     & ! [add] H.Yashiro 20110826
         input_size,        &
         layer_type,        &
         layername,         & ! [add] H.Yashiro 20110826
         nlayer,            & ! [add] H.Yashiro 20131030
         num_of_data,       &
         data_date,         &
         ddata_date,        & ! [Add] 12/02/01 T.Seiki
         opt_increment_date,& ! [Add] 12/02/01 T.Seiki
         fix_rec,           &
         opt_fix_rec,       &
         opt_monthly_cnst,  &
         opt_periodic_year, &
         defval

    real(RP) :: csec !!! [Add] T.Seiki, xxxxxx
    integer :: cdate(6)

    integer :: ierr
    integer :: np, im
    !---------------------------------------------------------------------------

    call calendar_ss2yh( cdate, ctime )

    !--- search the number of data files.
    max_extdata = 0
    rewind(ADM_CTL_FID)
    do
       read(ADM_CTL_FID, nml=nm_extdata, iostat=ierr)

       if ( ierr>0 ) then
          write(ADM_LOG_FID,*) 'Msg : Sub[extdata_setup]/Mod[mod_extdata]'
          write(ADM_LOG_FID,*) ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
          call ADM_proc_stop
       endif

       if(ierr < 0) exit

       max_extdata = max_extdata + 1
    enddo

    !--- allocation of data information array
    allocate( info(max_extdata) )

    ! -> [add&mod] H.Yashiro 20110826
    !--- read namelist, again.
    rewind(ADM_CTL_FID)
    do np=1, max_extdata

       !--- intialization
       fname             = ''
       dataname          = ''
       input_io_mode     = 'LEGACY' ! [add] H.Yashiro 20110826
       input_size        = 8
       layer_type        = ''
       layername         = 'ZSSFC1' ! [add] H.Yashiro 20110826
       num_of_data       = 1
       data_date(:,:)    = 1
       fix_rec           = 1
       opt_fix_rec       = .false.
       opt_monthly_cnst  = .false.
       opt_periodic_year = .false.
       defval            = 0.0_RP
       ! [Add] 12/02/01 T.Seiki
       ddata_date(:)     = 0
       opt_increment_date= .false.
       !
       !--- read namelist
       read(ADM_CTL_FID, nml=nm_extdata, iostat=ierr)

       if ( trim(layer_type) == 'ATM' ) then
          info(np)%kall = ADM_kall
       elseif( trim(layer_type) == 'SFC' ) then
          info(np)%kall = ADM_KNONE
       elseif( trim(layer_type) == 'NUM' ) then
          info(np)%kall = nlayer
       else
          write(ADM_LOG_FID,*) 'Msg : Sub[extdata_setup]/Mod[mod_extdata]'
          write(ADM_LOG_FID,*) 'xxx invlalid type of layer_type.', trim(layer_type)
          call ADM_proc_stop
       endif

       ! <- [add] H.Yashiro 20110826
       if ( input_io_mode == 'ADVANCED' ) then
          !--- Advanced I/O verifies;
          !---  info(np)%layername
          !---  info(np)%kall
          !
          !--- Advanced I/O overwrites;
          !---  info(np)%input_size
          !---  info(np)%num_of_data
          !---  info(np)%data_date
          !---  info(np)%data_rec(1)

          call FIO_seek( info(np)%data_rec(1), & !--- [out]
                         num_of_data,          & !--- [overwrite]
                         data_date,            & !--- [overwrite]
                         input_size,           & !--- [overwrite]
                         fname,                & !--- [in]
                         dataname,             & !--- [in]
                         layername,            & !--- [in]
                         1,                    & !--- [in]
                         info(np)%kall,        & !--- [in]
                         ctime,                & !--- [in]
                         cdate,                & !--- [in]
                         opt_periodic_year     ) !--- [in]

       elseif( input_io_mode == 'LEGACY' ) then
          ! [Add] 12/02/01 T.Seiki
          if (opt_increment_date)then
             do im=2, num_of_data
                data_date(1:6,im) = data_date(1:6,im-1) + ddata_date(1:6)
                call calendar_yh2ss( csec, data_date(:,im))
                call calendar_ss2yh( data_date(:,im), csec )
             enddo
          endif
          if (opt_periodic_year) then ! [add] C.Kodama 2010.07.26
             data_date(1,1:num_of_data) = cdate(1)
          endif
       else
          write(ADM_LOG_FID,*) 'xxx Invalid input_io_mode!', trim(input_io_mode)
          call ADM_proc_stop
       endif
       ! -> [add] H.Yashiro 20110826

       allocate( info(np)%data_date(6,num_of_data) )
       allocate( info(np)%data_time(  num_of_data) )
       allocate( info(np)%v   (ADM_gall,   info(np)%kall,ADM_lall   ,2) )
       allocate( info(np)%v_pl(ADM_gall_pl,info(np)%kall,ADM_lall_pl,2) )

       !--- store information
       info(np)%fname             = trim(fname)
       info(np)%dataname          = trim(dataname)
       info(np)%input_io_mode     = trim(input_io_mode) ! [add] H.Yashiro 20110826
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
          call ADM_proc_stop
       endif

       info(np)%data_date(:,1:num_of_data) = data_date(:,1:num_of_data)
       do im = 1, num_of_data
          call calendar_yh2ss( info(np)%data_time(im), info(np)%data_date(:,im) )
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

       elseif(info(np)%opt_periodic_year) then ! [add] C.Kodama 2010.07.26

          if( info(np)%input_io_mode == 'LEGACY' ) then
             info(np)%data_rec(1) = 1
             do im = 1, info(np)%num_of_data
                if ( ctime < info(np)%data_time(im) ) then
                   info(np)%data_rec(1) = im
                   exit
                endif
             enddo
          endif

          if( info(np)%data_rec(1) == 1 ) then
             info(np)%data_rec(2) = info(np)%num_of_data
          else
             info(np)%data_rec(2) = info(np)%data_rec(1)-1
          endif

       else !--- default

          if( info(np)%input_io_mode == 'LEGACY' ) then
             info(np)%data_rec(1) = 1
             do im = 1, info(np)%num_of_data
                if ( ctime < info(np)%data_time(im) ) then
                   info(np)%data_rec(1) = im
                   exit
                endif
             enddo
          endif

          if ( info(np)%data_rec(1) == 1 ) then
             write(ADM_LOG_FID,*) 'Msg : Sub[extdata_setup]/Mod[mod_extdata]'
             write(ADM_LOG_FID,*) '--- ERROR : data time is not consistent ', &
                       'with the simulation time! : ', trim(info(np)%dataname )
             call ADM_proc_stop
          else !--- default
             info(np)%data_rec(2) = info(np)%data_rec(1)-1
          endif

       endif

       !--- read initial data
       call data_read( np )

    enddo
    ! <- [add&mod] H.Yashiro 20110826

    !--- output information
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) 'Msg : Sub[extdata_setup]/Mod[mod_extdata]'
    write(ADM_LOG_FID,*) '===================================================='
    write(ADM_LOG_FID,*) '--- Number of maximum external data  : ',max_extdata
    do np = 1, max_extdata
       write(ADM_LOG_FID,'(1x,A,I4,A,A)') '--- variable [',np,'] : ', trim(info(np)%dataname)
    enddo
    write(ADM_LOG_FID,*) '===================================================='

    do np = 1, max_extdata
       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,'(1x,A,I4,A)') '============ External file NO. : ',np,' ============='
       write(ADM_LOG_FID,*) '--- fname             : ', trim(info(np)%fname)
       write(ADM_LOG_FID,*) '--- dataname          : ', trim(info(np)%dataname)
       write(ADM_LOG_FID,*) '--- input_io_mode     : ', trim(info(np)%input_io_mode) ! [add] H.Yashiro 20110826
       write(ADM_LOG_FID,*) '--- input_size        : ', info(np)%input_size
       write(ADM_LOG_FID,*) '--- layer_type        : ', info(np)%layer_type
       write(ADM_LOG_FID,*) '--- layername         : ', info(np)%layername ! [add] H.Yashiro 20110826
       write(ADM_LOG_FID,*) '--- num_of_data       : ', info(np)%num_of_data
       do im = 1, info(np)%num_of_data
          write(ADM_LOG_FID,'(1x,A,6(I4,1x))') '--- data_date         : ', info(np)%data_date(:,im)
       enddo
       write(ADM_LOG_FID,*) '--- opt_fix_rec       : ', info(np)%opt_fix_rec
       write(ADM_LOG_FID,*) '--- opt_monthly_cnst  : ', info(np)%opt_monthly_cnst
       write(ADM_LOG_FID,*) '--- opt_periodic_year : ', info(np)%opt_periodic_year
       write(ADM_LOG_FID,*) '--- defval            : ', info(np)%defval
       if ( info(np)%input_io_mode == 'ADVANCED' ) then
         write(ADM_LOG_FID,*) '--- first step        : ', info(np)%data_rec(1)
       endif
       write(ADM_LOG_FID,*) '====================================================='
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
    use mod_adm, only: &
      ADM_LOG_FID,     &
      ADM_gall_in,     &
      ADM_IopJop,      &
      ADM_GIoJo,       &
      ADM_proc_stop
    use mod_calendar, only: &
      calendar_yh2ss, &
      calendar_ss2yh
    implicit none

    real(RP),          intent(inout) :: gdata(:,:) ! data is inout to retain initilized value.
    character(len=*), intent(in)    :: DNAME      ! data name
    integer,          intent(in)    :: l_region
    real(RP),          intent(in)    :: ctime      ! current time
    logical,          intent(out)   :: eflag
    !
    real(RP) :: dt !-- delta t between two timestep data
    real(RP) :: wt !-- t weight of two timestep data

    !--- data ID
    integer :: np
    integer :: data_date_prev(6), cdate(6)
    real(RP) :: data_time_prev

    integer :: kall, gall
    integer :: im, n, k
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
       call ADM_proc_stop
    endif

    !--- update external data
    if (info(np)%opt_fix_rec) then
       !--- no time-advance
    else
       !--- data time advance
       if ( info(np)%opt_monthly_cnst ) then

          call calendar_ss2yh(cdate,ctime)

          if (cdate(2) /= info(np)%data_rec(1)) then
             info(np)%data_rec(1) = cdate(2)
             info(np)%data_rec(2) = info(np)%data_rec(1)

             call data_read( np )
          endif

       else !--- default

          if ( ctime > info(np)%data_time(info(np)%data_rec(1)) ) then
             !<-- current time pass the current data

             write(ADM_LOG_FID,*) '*** Update external data :',trim(info(np)%dataname)

             !--- increment of data_rec
             info(np)%data_rec(2) = info(np)%data_rec(1)
             info(np)%data_rec(1) = info(np)%data_rec(1) + 1
             if (       ( info(np)%data_rec(1) > info(np)%num_of_data ) &
                  .AND. ( .not. info(np)%opt_periodic_year )            ) then

                write(ADM_LOG_FID,*) 'xxx This run is over the land surface data range.'
                call ADM_proc_stop

             elseif( ( info(np)%data_rec(1) > info(np)%num_of_data ) .and. &
                     ( info(np)%opt_periodic_year ) ) then
                info(np)%data_rec(1) = info(np)%data_rec(1) - info(np)%num_of_data !-- 1
                info(np)%data_rec(2) = info(np)%num_of_data

                !--- re-setting of data_date and data_time
                info(np)%data_date(1,1:info(np)%num_of_data) = info(np)%data_date(1,1:info(np)%num_of_data) + 1
                do im = 1, info(np)%num_of_data
                   call calendar_yh2ss(info(np)%data_time(im),info(np)%data_date(:,im))
                enddo

                !--- output of message.
                write(ADM_LOG_FID,*) '*** data date is updated as follows.'
                do im = 1,info(np)%num_of_data
                   write(ADM_LOG_FID,*) ' ----- data_date(',im,') : ', info(np)%data_date(:,im)
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
       call calendar_yh2ss(data_time_prev,data_date_prev(:))

       dt = info(np)%data_time(info(np)%data_rec(1)) &
          - data_time_prev
       wt = ( ctime - data_time_prev ) / dt
    else !--- default
       dt = info(np)%data_time(info(np)%data_rec(1)) &
          - info(np)%data_time(info(np)%data_rec(2))
       wt = ( ctime - info(np)%data_time(info(np)%data_rec(2)) ) / dt
    endif

    do k = 1, kall
    do n = 1, ADM_gall_in
       gdata(n,k) = info(np)%v(ADM_IopJop(n,ADM_GIoJo),k,l_region,1) * (     wt) &
                  + info(np)%v(ADM_IopJop(n,ADM_GIoJo),k,l_region,2) * (1.0_RP-wt)
    enddo
    enddo

    return
  end subroutine extdata_update

  !-----------------------------------------------------------------------------
  subroutine data_read( np )
    use mod_adm, only : &
      ADM_LOG_FID
    use mod_fio, only : &
      FIO_input
    use mod_gtl, only : &
      GTL_input_var2_da
    use mod_comm, only : &
      COMM_var
    implicit none

    integer, intent(in) :: np

    integer :: n
    !---------------------------------------------------------------------------

    info(np)%v   (:,:,:,:) = info(np)%defval
    info(np)%v_pl(:,:,:,:) = info(np)%defval

    do n = 1, 2 !--- forward & backward
       ! <- [add] H.Yashiro 20110826
       if ( info(np)%input_io_mode == 'ADVANCED' ) then

          call FIO_input( info(np)%v(:,:,:,n), &
                          info(np)%fname,      &
                          info(np)%dataname,   &
                          info(np)%layername,  &
                          1,info(np)%kall,     &
                          info(np)%data_rec(n) )

       elseif( info(np)%input_io_mode == 'LEGACY' ) then
       ! -> [add] H.Yashiro 20110826

          call GTL_input_var2_da( trim(info(np)%fname),            &
                                  info(np)%v(:,:,:,n),             &
                                  1, info(np)%kall,                &
                                  recnum=info(np)%data_rec(n),     &
                                  input_size = info(np)%input_size )

       endif
    enddo

    call COMM_var( info(np)%v, info(np)%v_pl, info(np)%kall, 2 )

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '===== READ EXTERNAL DATA ============================'
    write(ADM_LOG_FID,'(1x,A,I4,A,A)') '--- variable [',np,'] : ', trim(info(np)%dataname)
    write(ADM_LOG_FID,*) '--- forward  data step(record) number : ',info(np)%data_rec(1)
    write(ADM_LOG_FID,*) '--- backward data step(record) number : ',info(np)%data_rec(2)
    write(ADM_LOG_FID,*) '====================================================='

    return
  end subroutine data_read
  !-----------------------------------------------------------------------------

end module mod_extdata
!-------------------------------------------------------------------------------
