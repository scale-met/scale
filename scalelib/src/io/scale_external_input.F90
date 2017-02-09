!-------------------------------------------------------------------------------
!> module EXTERNAL INPUT
!!
!! @par Description
!!          External file input module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_external_input
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_land_grid_index
  use scale_urban_grid_index
  use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: EXTIN_setup
  public :: EXTIN_regist
  public :: EXTIN_update

  interface EXTIN_update
     module procedure EXTIN_update_1D
     module procedure EXTIN_update_2D
     module procedure EXTIN_update_3D
  end interface EXTIN_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: EXTIN_time_advance

  abstract interface
     subroutine get_dims_1D( &
          dim1_max, &
          dim1_S,   &
          dim1_E,   &
          varname,  &
          axistype  )
       integer,          intent(out) :: dim1_max
       integer,          intent(out) :: dim1_S
       integer,          intent(out) :: dim1_E
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: axistype     ! axis type (Z/X/Y)
     end subroutine get_dims_1D
     subroutine get_dims_2D( &
          dim1_max,  &
          dim1_S,    &
          dim1_E,    &
          dim2_max,  &
          dim2_S,    &
          dim2_E,    &
          transpose, &
          varname,   &
          axistype   )
       integer,          intent(out) :: dim1_max
       integer,          intent(out) :: dim1_S
       integer,          intent(out) :: dim1_E
       integer,          intent(out) :: dim2_max
       integer,          intent(out) :: dim2_S
       integer,          intent(out) :: dim2_E
       logical,          intent(out) :: transpose
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: axistype     ! axis type (XY/XZ/ZX)
     end subroutine get_dims_2D
     subroutine get_dims_3D( &
          dim1_max,  &
          dim1_S,    &
          dim1_E,    &
          dim2_max,  &
          dim2_S,    &
          dim2_E,    &
          dim3_max,  &
          dim3_S,    &
          dim3_E,    &
          transpose, &
          varname,   &
          axistype   )
       integer,          intent(out) :: dim1_max
       integer,          intent(out) :: dim1_S
       integer,          intent(out) :: dim1_E
       integer,          intent(out) :: dim2_max
       integer,          intent(out) :: dim2_S
       integer,          intent(out) :: dim2_E
       integer,          intent(out) :: dim3_max
       integer,          intent(out) :: dim3_S
       integer,          intent(out) :: dim3_E
       logical,          intent(out) :: transpose
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: axistype     ! axis type (ZXY/XYZ/Land/Urban)
     end subroutine get_dims_3D
  end interface
  procedure(get_dims_1D), pointer :: EXTIN_get_dims_1D => NULL()
  procedure(get_dims_2D), pointer :: EXTIN_get_dims_2D => NULL()
  procedure(get_dims_3D), pointer :: EXTIN_get_dims_3D => NULL()
  private :: EXTIN_get_dims_1D
  private :: EXTIN_get_dims_2D
  private :: EXTIN_get_dims_3D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: I_prev = 1 !< [index] previous
  integer, private, parameter :: I_next = 2 !< [index] next

  integer, private, parameter :: I_periodic_year  = 1
  integer, private, parameter :: I_periodic_month = 2
  integer, private, parameter :: I_periodic_day   = 3

  integer, private, parameter :: EXTIN_item_limit = 1000  !< limit of item
  integer, private, parameter :: EXTIN_step_limit = 10000 !< limit of steps          for each item
  integer, private, parameter :: EXTIN_dim_limit  = 3     !< limit of dimension rank for each item

  type, private :: EXTIN_itemcontainer
     character(len=H_SHORT) :: varname                   !< variable name
     integer                :: fid                       !< file id
     integer                :: ndim                      !< number of dimensions
     integer                :: dim_size(EXTIN_dim_limit) !< size of dimension (z,x,y)
     integer, allocatable   :: dim_start(:)              !< start index
     integer                :: step_num                  !< size of time dimension
     real(DP), allocatable  :: time(:)                   !< time of each step [sec]
     logical                :: fixed_step                !< fix step position?
     integer                :: flag_periodic             !< treat as periodic data? (0:no 1:yearly 2:monthly 3:daily)
     real(RP)               :: offset                    !< offset value (default is set to 0)
     integer                :: data_steppos(2)           !< current step position to read (1:previous,2:next)
     real(RP), allocatable  :: value(:,:,:,:)            !< data value                    (1:previous,2:next)
     character(len=H_SHORT) :: axistype                  ! axis type
     logical                :: transpose                 !< true: xyz, false: zxy
  end type EXTIN_itemcontainer


  integer,                   private :: EXTIN_item_count = 0     !< number of item to output
  type(EXTIN_itemcontainer), private :: EXTIN_item(EXTIN_item_limit)            !< item to output

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine EXTIN_setup( model )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_external_input_rm, only: &
       EXTIN_RM_get_dims_1D, &
       EXTIN_RM_get_dims_2D, &
       EXTIN_RM_get_dims_3D
    implicit none
    character(len=*), intent(in) :: model

    character(len=H_LONG)  :: basename
    character(len=H_SHORT) :: varname
    character(len=H_SHORT) :: axistype
    integer                :: step_limit            ! limit number for reading data
    integer                :: step_fixed            ! fixed step position to read
    logical                :: enable_periodic_year  ! treat as yearly               periodic data?
    logical                :: enable_periodic_month ! treat as yearly,monthly       periodic data?
    logical                :: enable_periodic_day   ! treat as yearly,monthly,daily periodic data?
    real(RP)               :: offset
    real(RP)               :: defval
    logical                :: check_coordinates

    namelist /EXTITEM/ &
       basename,              &
       varname,               &
       axistype,              &
       step_limit,            &
       step_fixed,            &
       enable_periodic_year,  &
       enable_periodic_month, &
       enable_periodic_day,   &
       offset,                &
       defval,                &
       check_coordinates

    integer  :: count
    integer  :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[EXTIN] / Categ[ATMOS-RM IO] / Origin[SCALElib]'

    select case ( model )
    case ( 'RM' )
       EXTIN_get_dims_1D => EXTIN_RM_get_dims_1D
       EXTIN_get_dims_2D => EXTIN_RM_get_dims_2D
       EXTIN_get_dims_3D => EXTIN_RM_get_dims_3D
!    case ( 'GM' )
    case default
       write(*,*) 'xxx EXTIN is not support for the model: ', trim(model)
       call PRC_MPIstop
    end select

    ! count external data from namelist
    rewind(IO_FID_CONF)
    do count = 1, EXTIN_item_limit

       ! set default
       step_limit            = EXTIN_step_limit
       basename              = ''
       varname               = ''
       axistype              = ''
       step_fixed            = -1
       enable_periodic_year  = .false.
       enable_periodic_month = .false.
       enable_periodic_day   = .false.
       offset                = 0.0_RP
       defval                = UNDEF
       check_coordinates     = .false.

       ! read namelist
       read(IO_FID_CONF,nml=EXTITEM,iostat=ierr)

       if( ierr < 0 ) then !--- no more items
          exit
       elseif( ierr > 0 ) then !--- fatal error
          write(*,*) 'xxx Not appropriate names in namelist EXTITEM. Check!', count
          call PRC_MPIstop
       endif

       if( IO_NML .AND. IO_FID_NML /= IO_FID_LOG ) write(IO_FID_NML,nml=EXTITEM)

       call EXTIN_regist(          &
            basename,              &
            varname,               &
            axistype,              &
            enable_periodic_year,  &
            enable_periodic_month, &
            enable_periodic_day,   &
            step_fixed,            &
            offset,                &
            defval,                &
            check_coordinates,     &
            step_limit             )

    enddo

    return
  end subroutine EXTIN_setup

  !-----------------------------------------------------------------------------
  !> Regist data
  subroutine EXTIN_regist( &
       basename,              &
       varname,               &
       axistype,              &
       enable_periodic_year,  &
       enable_periodic_month, &
       enable_periodic_day,   &
       step_fixed,            &
       offset,                &
       defval,                &
       check_coordinates,     &
       step_limit             )
    use gtool_file_h, only: &
       File_FREAD
    use gtool_file, only: &
       FileOpen, &
       FileGetAllDataInfo, &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_fileio, only: &
       FILEIO_check_coordinates
    use scale_calendar, only: &
       CALENDAR_adjust_daysec,  &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_combine_daysec, &
       CALENDAR_CFunits2sec,    &
       I_year,                  &
       I_month,                 &
       I_day
    use scale_time, only: &
       TIME_STARTDAYSEC, &
       TIME_NOWDAYSEC, &
       TIME_OFFSET_year
    implicit none

    character(len=*) , intent(in) :: basename
    character(len=*) , intent(in) :: varname
    character(len=*) , intent(in) :: axistype
    integer          , intent(in) :: step_fixed            ! fixed step position to read
    logical          , intent(in) :: enable_periodic_year  ! treat as yearly               periodic data?
    logical          , intent(in) :: enable_periodic_month ! treat as yearly,monthly       periodic data?
    logical          , intent(in) :: enable_periodic_day   ! treat as yearly,monthly,daily periodic data?
    real(RP)         , intent(in) :: offset
    real(RP)         , intent(in) :: defval
    logical, optional, intent(in) :: check_coordinates
    integer, optional, intent(in) :: step_limit            ! limit number for reading data

    integer                :: step_nmax
    character(len=H_LONG)  :: description
    character(len=H_SHORT) :: unit
    integer                :: datatype
    integer                :: dim_rank
    character(len=H_SHORT) :: dim_name  (3)
    integer                :: dim_size  (3)
    real(DP)               :: time_start(EXTIN_step_limit)
    real(DP)               :: time_end  (EXTIN_step_limit)
    character(len=H_MID)   :: time_units

    integer  :: datadate(6)   !< date
    real(DP) :: datasubsec    !< subsecond
    integer  :: dataday       !< absolute day
    real(DP) :: datasec       !< absolute second
    integer  :: offset_year   !< offset year

    integer  :: dim1_max, dim1_S, dim1_E
    integer  :: dim2_max, dim2_S, dim2_E
    integer  :: dim3_max, dim3_S, dim3_E

    integer  :: step_limit_

    integer  :: fid
    integer  :: nid, n

    if ( present(step_limit) ) then
       if ( step_limit > 0 ) then
          step_limit_ = step_limit
       else
          step_limit_ = EXTIN_step_limit
       end if
    else
       step_limit_ = EXTIN_step_limit
    end if

    do nid = 1, EXTIN_item_count
       if ( EXTIN_item(nid)%varname  == varname ) then
          write(*,*) 'xxx Data is already registered! basename,varname = ', trim(basename), ', ', trim(varname)
          call PRC_MPIstop
       end if
    end do

    EXTIN_item_count = EXTIN_item_count + 1

    if ( EXTIN_item_count > EXTIN_item_limit ) then
       write(*,*) 'xxx Number of EXT data exceedes the limit', EXTIN_item_count, EXTIN_item_limit
       call PRC_MPIstop
    end if


    call FileOpen( fid, basename, File_FREAD, myrank=PRC_myrank )

    ! read from file
    call FileGetAllDatainfo( step_limit_,               & ! [IN]
                             EXTIN_dim_limit,           & ! [IN]
                             fid,                       & ! [IN]
                             varname,                   & ! [IN]
                             step_nmax,                 & ! [OUT]
                             description,               & ! [OUT]
                             unit,                      & ! [OUT]
                             datatype,                  & ! [OUT]
                             dim_rank,                  & ! [OUT]
                             dim_name  (:),             & ! [OUT]
                             dim_size  (:),             & ! [OUT]
                             time_start(1:step_limit_), & ! [OUT]
                             time_end  (1:step_limit_), & ! [OUT]
                             time_units                 ) ! [OUT]

    if ( step_nmax == 0 ) then
       write(*,*) 'xxx Data not found! basename,varname = ', trim(basename), ', ', trim(varname)
       call PRC_MPIstop
    endif

    do n = dim_rank+1, 3
       dim_size(n) = 1
    end do

    nid = EXTIN_item_count

    ! setup item
    EXTIN_item(nid)%fid         = fid
    EXTIN_item(nid)%varname     = varname
    EXTIN_item(nid)%dim_size(:) = dim_size(:)
    EXTIN_item(nid)%step_num    = step_nmax

    if ( enable_periodic_day ) then
       EXTIN_item(nid)%flag_periodic = I_periodic_day
    elseif( enable_periodic_month ) then
       EXTIN_item(nid)%flag_periodic = I_periodic_month
    elseif( enable_periodic_year ) then
       EXTIN_item(nid)%flag_periodic = I_periodic_year
    else
       EXTIN_item(nid)%flag_periodic = 0
    endif

    allocate( EXTIN_item(nid)%value(dim_size(1),dim_size(2),dim_size(3),2) )
    EXTIN_item(nid)%value(:,:,:,:) = defval
    EXTIN_item(nid)%offset         = offset

    allocate( EXTIN_item(nid)%time(step_nmax) )
    EXTIN_item(nid)%time(:) = 0.0_DP

    do n = 1, EXTIN_item(nid)%step_num
       EXTIN_item(nid)%time(n) = CALENDAR_CFunits2sec( time_end(n), time_units, TIME_OFFSET_year, TIME_STARTDAYSEC )
    enddo

    if ( EXTIN_item(nid)%step_num == 1 ) then

       EXTIN_item(nid)%fixed_step           = .true.
       EXTIN_item(nid)%data_steppos(I_prev) = 1
       EXTIN_item(nid)%data_steppos(I_next) = 1

    else if ( step_fixed > 0 ) then ! fixed time step mode

       EXTIN_item(nid)%fixed_step           = .true.
       EXTIN_item(nid)%data_steppos(I_prev) = step_fixed
       EXTIN_item(nid)%data_steppos(I_next) = step_fixed

    else

       EXTIN_item(nid)%fixed_step = .false.

       ! seek start position
       EXTIN_item(nid)%data_steppos(I_next) = 1
       do n = 1, EXTIN_item(nid)%step_num
          if ( EXTIN_item(nid)%time(n) > TIME_NOWDAYSEC ) exit
          EXTIN_item(nid)%data_steppos(I_next) = n + 1
       enddo

       EXTIN_item(nid)%data_steppos(I_prev) = EXTIN_item(nid)%data_steppos(I_next) - 1

       if ( EXTIN_item(nid)%flag_periodic > 0 ) then ! periodic time step mode

          if ( EXTIN_item(nid)%data_steppos(I_next) == 1 ) then ! between first-1 and first

             ! first-1 = last
             EXTIN_item(nid)%data_steppos(I_prev) = EXTIN_item(nid)%step_num

          elseif( EXTIN_item(nid)%data_steppos(I_next) == EXTIN_item(nid)%step_num+1 ) then ! between last and last+1

             ! last+1 = first
             EXTIN_item(nid)%data_steppos(I_next) = 1

             ! update data time in periodic condition
             do n = 1, EXTIN_item(nid)%step_num
                dataday     = 0
                datasec     = EXTIN_item(nid)%time(n)
                offset_year = 0
                call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

                call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                           datasubsec,  & ! [OUT]
                                           dataday,     & ! [IN]
                                           datasec,     & ! [IN]
                                           offset_year  ) ! [IN]

                if    ( EXTIN_item(nid)%flag_periodic == I_periodic_day   ) then
                   datadate(I_day)   = datadate(I_day)   + 1
                elseif( EXTIN_item(nid)%flag_periodic == I_periodic_month ) then
                   datadate(I_month) = datadate(I_month) + 1
                elseif( EXTIN_item(nid)%flag_periodic == I_periodic_year  ) then
                   datadate(I_year)  = datadate(I_year)  + 1
                endif

                call CALENDAR_date2daysec( dataday,     & ! [OUT]
                                           datasec,     & ! [OUT]
                                           datadate(:), & ! [IN]
                                           datasubsec,  & ! [IN]
                                           offset_year  ) ! [IN]

                EXTIN_item(nid)%time(n) = CALENDAR_combine_daysec( dataday, datasec )
             enddo

             if( IO_L ) write(IO_FID_LOG,*) '*** data time is updated.'
          endif

       else ! normal mode

          if (      EXTIN_item(nid)%data_steppos(I_next) == 1                          &
               .OR. EXTIN_item(nid)%data_steppos(I_next) == EXTIN_item(nid)%step_num+1 ) then
             write(*,*) 'xxx Current time is out of period of external data! ', trim(varname)
             call PRC_MPIstop
          endif

       endif

    endif

    !--- read first data
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Initial read of external data : ', trim(varname)

    if (       dim_size(1) >= 1 &
         .AND. dim_size(2) == 1 &
         .AND. dim_size(3) == 1 ) then ! 1D

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 1D var                   : ', trim(varname)

       call EXTIN_get_dims_1D( &
            dim1_max, dim1_S, dim1_E, & ! (out)
            varname, axistype         ) ! (in)

       EXTIN_item(nid)%ndim      = 1
       EXTIN_item(nid)%transpose = .false.
       allocate( EXTIN_item(nid)%dim_start(1) )
       EXTIN_item(nid)%dim_start(1) = dim1_S

       if ( dim1_max /= dim_size(1) ) then
          write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
          write(*,*) 'xxx dim 1 (data,requested)    : ', dim_size(1), dim1_max
          call PRC_MPIstop
       endif

       ! read prev
       call FileRead( EXTIN_item(nid)%value(:,1,1,I_prev), &
                      EXTIN_item(nid)%fid,                 &
                      EXTIN_item(nid)%varname,             &
                      EXTIN_item(nid)%data_steppos(I_prev) )
       ! read next
       call FileRead( EXTIN_item(nid)%value(:,1,1,I_next), &
                      EXTIN_item(nid)%fid,                 &
                      EXTIN_item(nid)%varname,             &
                      EXTIN_item(nid)%data_steppos(I_next) )

    elseif (       dim_size(1) >= 1 &
             .AND. dim_size(2) >  1 &
             .AND. dim_size(3) == 1 ) then ! 2D

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var                   : ', trim(varname)

       call EXTIN_get_dims_2D( &
            dim1_max, dim1_S, dim1_E,  & ! (out)
            dim2_max, dim2_S, dim2_E,  & ! (out)
            EXTIN_item(nid)%transpose, & ! (out)
            varname, axistype          ) ! (in)

       EXTIN_item(nid)%ndim      = 2
       allocate( EXTIN_item(nid)%dim_start(2) )
       EXTIN_item(nid)%dim_start(1) = dim1_S
       EXTIN_item(nid)%dim_start(2) = dim2_S

       if (      dim1_max /= dim_size(1) &
            .OR. dim2_max /= dim_size(2) ) then
          write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
          write(*,*) 'xxx dim 1 (data,requested)    : ', dim_size(1), dim1_max
          write(*,*) 'xxx dim 2 (data,requested)    : ', dim_size(2), dim2_max
          call PRC_MPIstop
       endif

       ! read prev
       call FileRead( EXTIN_item(nid)%value(:,:,1,I_prev), &
                      EXTIN_item(nid)%fid,                 &
                      EXTIN_item(nid)%varname,             &
                      EXTIN_item(nid)%data_steppos(I_prev) )
       ! read next
       call FileRead( EXTIN_item(nid)%value(:,:,1,I_next), &
                      EXTIN_item(nid)%fid,                 &
                      EXTIN_item(nid)%varname,             &
                      EXTIN_item(nid)%data_steppos(I_next) )

    elseif (       dim_size(1) >= 1 &
             .AND. dim_size(2) >  1 &
             .AND. dim_size(3) >  1 ) then ! 3D

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var                   : ', trim(varname)

       call EXTIN_get_dims_3D( &
            dim1_max, dim1_S, dim1_E,  & ! (out)
            dim2_max, dim2_S, dim2_E,  & ! (out)
            dim3_max, dim3_S, dim3_E,  & ! (out)
            EXTIN_item(nid)%transpose, & ! (out)
            varname, axistype          ) ! (in)

       EXTIN_item(nid)%ndim      = 3
       allocate( EXTIN_item(nid)%dim_start(3) )
       EXTIN_item(nid)%dim_start(1) = dim1_S
       EXTIN_item(nid)%dim_start(2) = dim2_S
       EXTIN_item(nid)%dim_start(3) = dim3_S

       if (      dim1_max /= dim_size(1) &
            .OR. dim2_max /= dim_size(2) &
            .OR. dim3_max /= dim_size(3) ) then
          write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
          write(*,*) 'xxx dim 1 (data,requested)    : ', dim_size(1), dim1_max
          write(*,*) 'xxx dim 2 (data,requested)    : ', dim_size(2), dim2_max
          write(*,*) 'xxx dim 3 (data,requested)    : ', dim_size(3), dim3_max
          call PRC_MPIstop
       endif

       ! read prev
       call FileRead( EXTIN_item(nid)%value(:,:,:,I_prev), &
                      EXTIN_item(nid)%fid,                 &
                      EXTIN_item(nid)%varname,             &
                      EXTIN_item(nid)%data_steppos(I_prev) )
       ! read next
       call FileRead( EXTIN_item(nid)%value(:,:,:,I_next), &
                      EXTIN_item(nid)%fid,                 &
                      EXTIN_item(nid)%varname,             &
                      EXTIN_item(nid)%data_steppos(I_next) )

    else
       write(*,*) 'xxx Unexpected dimsize: ', dim_size(:)
       call PRC_MPIstop
    endif

    if ( present(check_coordinates) ) then
       if ( check_coordinates ) &
            call FILEIO_check_coordinates( fid, &
                    atmos     = EXTIN_item(nid)%ndim==3, &
                    transpose = EXTIN_item(nid)%transpose )
    end if

    return
  end subroutine EXTIN_regist

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine EXTIN_update_1D( &
       var,          &
       varname,      &
       current_time, &
       error         )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:)       ! variable
    character(len=*), intent(in)  :: varname      ! item name
    real(DP),         intent(in)  :: current_time ! current time
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile

    integer  :: n
    integer  :: n1
    integer  :: nn1
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, EXTIN_item_count
       if( varname == EXTIN_item(n)%varname ) nid = n
    enddo

    if ( nid == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Variable was not registered: ', trim(varname)
       return
    end if

    if ( EXTIN_item(nid)%ndim /= 1 ) then
       write(*,*) 'xxx Data is not 1D var: ', trim(EXTIN_item(nid)%varname)
       call PRC_MPIstop
    end if

    call EXTIN_time_advance( nid,          & ! [IN]
                             current_time, & ! [IN]
                             weight,       & ! [OUT]
                             do_readfile   ) ! [OUT]

    if ( do_readfile ) then
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 1D var           : ', trim(EXTIN_item(nid)%varname)

       ! next -> prev
       EXTIN_item(nid)%value(:,:,:,I_prev) = EXTIN_item(nid)%value(:,:,:,I_next)

       ! read next
       call FileRead( EXTIN_item(nid)%value(:,1,1,I_next), & ! [OUT]
                      EXTIN_item(nid)%fid,                 & ! [IN]
                      EXTIN_item(nid)%varname,             & ! [IN]
                      EXTIN_item(nid)%data_steppos(I_next) )
    endif

    ! store data with weight
    do n1 = 1, EXTIN_item(nid)%dim_size(1)
    nn1 = n1 + EXTIN_item(nid)%dim_start(1) - 1

       var(nn1) = ( 1.0_RP-weight ) * EXTIN_item(nid)%value(n1,1,1,I_prev) &
                + (        weight ) * EXTIN_item(nid)%value(n1,1,1,I_next)
    enddo

    error = .false.

    return
  end subroutine EXTIN_update_1D

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine EXTIN_update_2D( &
       var,          &
       varname,      &
       current_time, &
       error         )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:)     ! variable
    character(len=*), intent(in)  :: varname      ! item name
    real(DP),         intent(in)  :: current_time ! current time
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile

    integer  :: n
    integer  :: n1, n2
    integer  :: nn1, nn2
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, EXTIN_item_count
       if( varname == EXTIN_item(n)%varname ) nid = n
    enddo

    if ( nid == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx variable was not registered: ', trim(varname)
       return
    end if

    call EXTIN_time_advance( nid,          & ! [IN]
                             current_time, & ! [IN]
                             weight,       & ! [OUT]
                             do_readfile   ) ! [OUT]

    if ( do_readfile ) then

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var           : ', trim(EXTIN_item(nid)%varname)
       ! next -> prev
       EXTIN_item(nid)%value(:,:,:,I_prev) = EXTIN_item(nid)%value(:,:,:,I_next)

       ! read next
       call FileRead( EXTIN_item(nid)%value(:,:,1,I_next), & ! [OUT]
                      EXTIN_item(nid)%fid,                 & ! [IN]
                      EXTIN_item(nid)%varname,             & ! [IN]
                      EXTIN_item(nid)%data_steppos(I_next) ) ! [IN]
    endif

    if ( EXTIN_item(nid)%transpose ) then
       ! store data with weight (x,z)->(z,x)
       do n1 = 1, EXTIN_item(nid)%dim_size(1)
       nn1 = n1 + EXTIN_item(nid)%dim_start(1) - 1
       do n2 = 1, EXTIN_item(nid)%dim_size(2)

       nn2 = n2 + EXTIN_item(nid)%dim_start(2) - 1
          var(nn2,nn1) = ( 1.0_RP-weight ) * EXTIN_item(nid)%value(n1,n2,1,I_prev) &
                       + (        weight ) * EXTIN_item(nid)%value(n1,n2,1,I_next)
       enddo
       enddo
    else
       ! store data with weight
       do n2 = 1, EXTIN_item(nid)%dim_size(2)
       nn2 = n2 + EXTIN_item(nid)%dim_start(2) - 1
       do n1 = 1, EXTIN_item(nid)%dim_size(1)
       nn1 = n1 + EXTIN_item(nid)%dim_start(1) - 1

          var(nn1,nn2) = ( 1.0_RP-weight ) * EXTIN_item(nid)%value(n1,n2,1,I_prev) &
                       + (        weight ) * EXTIN_item(nid)%value(n1,n2,1,I_next)
       enddo
       enddo
    end if

    error = .false.

    return
  end subroutine EXTIN_update_2D

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine EXTIN_update_3D( &
       var,          &
       varname,      &
       current_time, &
       error         )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:,:)   ! variable
    character(len=*), intent(in)  :: varname      ! item name
    real(DP),         intent(in)  :: current_time ! current time
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile

    integer  :: n
    integer  :: n1, n2, n3
    integer  :: nn1, nn2, nn3

    integer :: tmp
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, EXTIN_item_count
       if( varname == EXTIN_item(n)%varname ) nid = n
    enddo

    if ( nid == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx variable was not registered: ', trim(varname)
       return
    end if

    call EXTIN_time_advance( nid,          & ! [IN]
                             current_time, & ! [IN]
                             weight,       & ! [OUT]
                             do_readfile   ) ! [OUT]

    if ( do_readfile ) then
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var           : ', trim(EXTIN_item(nid)%varname)
       ! next -> prev
       EXTIN_item(nid)%value(:,:,:,I_prev) = EXTIN_item(nid)%value(:,:,:,I_next)

       ! read next
       call FileRead( EXTIN_item(nid)%value(:,:,:,I_next), & ! [OUT]
                      EXTIN_item(nid)%fid,                 & ! [IN]
                      EXTIN_item(nid)%varname,             & ! [IN]
                      EXTIN_item(nid)%data_steppos(I_next) ) ! [IN]
    endif

    if ( EXTIN_item(nid)%transpose ) then
       ! store data with weight (x,y,z)->(z,x,y)
       do n2 = 1, EXTIN_item(nid)%dim_size(2)
       nn2 = n2 + EXTIN_item(nid)%dim_start(2) - 1
       do n1 = 1, EXTIN_item(nid)%dim_size(1)
       nn1 = n1 + EXTIN_item(nid)%dim_start(1) - 1
       do n3 = 1, EXTIN_item(nid)%dim_size(3)
       nn3 = n3 + EXTIN_item(nid)%dim_start(3) - 1

          var(nn3,nn1,nn2) = ( 1.0_RP-weight ) * EXTIN_item(nid)%value(n1,n2,n3,I_prev) &
                           + (        weight ) * EXTIN_item(nid)%value(n1,n2,n3,I_next)
       enddo
       enddo
       enddo
    else
       ! store data with weight (z,x,y)->(z,x,y)
       do n3 = 1, EXTIN_item(nid)%dim_size(3)
       nn3 = n3 + EXTIN_item(nid)%dim_start(3) - 1
       do n2 = 1, EXTIN_item(nid)%dim_size(2)
       nn2 = n2 + EXTIN_item(nid)%dim_start(2) - 1
       do n1 = 1, EXTIN_item(nid)%dim_size(1)
       nn1 = n1 + EXTIN_item(nid)%dim_start(1) - 1

          var(nn1,nn2,nn3) = ( 1.0_RP-weight ) * EXTIN_item(nid)%value(n1,n2,n3,I_prev) &
                           + (        weight ) * EXTIN_item(nid)%value(n1,n2,n3,I_next)
       enddo
       enddo
       enddo
    end if

    error = .false.

    return
  end subroutine EXTIN_update_3D

  !-----------------------------------------------------------------------------
  !> Check time to read and calc weight
  subroutine EXTIN_time_advance( &
       nid,          &
       current_time, &
       weight,       &
       do_readfile   )
    use scale_calendar, only: &
       CALENDAR_adjust_daysec,  &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_combine_daysec, &
       I_year,                  &
       I_month,                 &
       I_day
    implicit none

    integer,  intent(in)  :: nid          ! item id
    real(DP), intent(in)  :: current_time ! current time
    real(RP), intent(out) :: weight       ! weight
    logical,  intent(out) :: do_readfile  ! read new data at this time?

    integer  :: datadate(6)   !< date
    real(DP) :: datasubsec    !< subsecond
    integer  :: dataday       !< absolute day
    real(DP) :: datasec       !< absolute second
    integer  :: offset_year   !< offset year

    real(DP) :: prev_time, next_time
    integer  :: t
    !---------------------------------------------------------------------------

    do_readfile = .false.

    if ( EXTIN_item(nid)%fixed_step ) then
       !--- no time-advance
    else
       ! time is passed?
       if ( current_time > EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_next) ) ) then

          do_readfile = .true.

          if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Update external input : ', trim(EXTIN_item(nid)%varname)

          ! update step position
          EXTIN_item(nid)%data_steppos(I_prev) = EXTIN_item(nid)%data_steppos(I_prev) + 1
          EXTIN_item(nid)%data_steppos(I_next) = EXTIN_item(nid)%data_steppos(I_next) + 1

          if ( EXTIN_item(nid)%flag_periodic > 0 ) then ! periodic time step mode

             if ( EXTIN_item(nid)%data_steppos(I_next) == EXTIN_item(nid)%step_num+1 ) then

                ! last+1 = first
                EXTIN_item(nid)%data_steppos(I_next) = 1

                ! update data time in periodic condition
                do t = 1, EXTIN_item(nid)%step_num
                   dataday     = 0
                   datasec     = EXTIN_item(nid)%time(t)
                   offset_year = 0
                   call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

                   call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                              datasubsec,  & ! [OUT]
                                              dataday,     & ! [IN]
                                              datasec,     & ! [IN]
                                              offset_year  ) ! [IN]

                   if    ( EXTIN_item(nid)%flag_periodic == I_periodic_day   ) then
                      datadate(I_day)   = datadate(I_day)   + 1
                   elseif( EXTIN_item(nid)%flag_periodic == I_periodic_month ) then
                      datadate(I_month) = datadate(I_month) + 1
                   elseif( EXTIN_item(nid)%flag_periodic == I_periodic_year  ) then
                      datadate(I_year)  = datadate(I_year)  + 1
                   endif

                   call CALENDAR_date2daysec( dataday,     & ! [IN]
                                              datasec,     & ! [IN]
                                              datadate(:), & ! [OUT]
                                              datasubsec,  & ! [OUT]
                                              offset_year  ) ! [IN]

                   EXTIN_item(nid)%time(t) = CALENDAR_combine_daysec( dataday, datasec )
                enddo
             endif

          else ! normal mode

             if ( EXTIN_item(nid)%data_steppos(I_next) == EXTIN_item(nid)%step_num+1 ) then
                write(*,*) 'xxx Current time is out of period of external data! '
                call PRC_MPIstop
             endif

          endif

       endif

    endif

    ! calc weight
    if ( EXTIN_item(nid)%fixed_step ) then

       weight = 0.0_RP

    elseif( EXTIN_item(nid)%data_steppos(I_next) == 1 ) then !<--- this case is only periodic one.

       dataday     = 0
       datasec     = EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_prev) )
       offset_year = 0
       call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

       call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                  datasubsec,  & ! [OUT]
                                  dataday,     & ! [IN]
                                  datasec,     & ! [IN]
                                  offset_year  ) ! [IN]

       if    ( EXTIN_item(nid)%flag_periodic == I_periodic_day   ) then
          datadate(I_day)   = datadate(I_day)   - 1
       elseif( EXTIN_item(nid)%flag_periodic == I_periodic_month ) then
          datadate(I_month) = datadate(I_month) - 1
       elseif( EXTIN_item(nid)%flag_periodic == I_periodic_year  ) then
          datadate(I_year)  = datadate(I_year)  - 1
       endif

       call CALENDAR_date2daysec( dataday,     & ! [IN]
                                  datasec,     & ! [IN]
                                  datadate(:), & ! [OUT]
                                  datasubsec,  & ! [OUT]
                                  offset_year  ) ! [IN]

       prev_time = CALENDAR_combine_daysec( dataday, datasec )

       next_time = EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_next) )

       weight = ( current_time - prev_time ) &
              / ( next_time    - prev_time )

    else

       prev_time = EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_prev) )
       next_time = EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_next) )

       weight = ( current_time - prev_time ) &
              / ( next_time    - prev_time )

    endif

    return
  end subroutine EXTIN_time_advance

end module scale_external_input
