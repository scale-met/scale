!-------------------------------------------------------------------------------
!> module file / external_input
!!
!! @par Description
!!          External file input module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_file_external_input
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
  public :: FILE_EXTERNAL_INPUT_setup
  public :: FILE_EXTERNAL_INPUT_regist
  public :: FILE_EXTERNAL_INPUT_update

  interface FILE_EXTERNAL_INPUT_update
     module procedure FILE_EXTERNAL_INPUT_update_1D
     module procedure FILE_EXTERNAL_INPUT_update_2D
     module procedure FILE_EXTERNAL_INPUT_update_3D
  end interface FILE_EXTERNAL_INPUT_update

  abstract interface
     subroutine get_dims1D( &
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
     end subroutine get_dims1D

     subroutine get_dims2D( &
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
     end subroutine get_dims2D

     subroutine get_dims3D( &
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
     end subroutine get_dims3D
  end interface

  procedure(get_dims1D), pointer :: FILE_EXTERNAL_INPUT_get_dims1D => NULL()
  procedure(get_dims2D), pointer :: FILE_EXTERNAL_INPUT_get_dims2D => NULL()
  procedure(get_dims3D), pointer :: FILE_EXTERNAL_INPUT_get_dims3D => NULL()
  public :: FILE_EXTERNAL_INPUT_get_dims1D
  public :: FILE_EXTERNAL_INPUT_get_dims2D
  public :: FILE_EXTERNAL_INPUT_get_dims3D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: FILE_EXTERNAL_INPUT_file_limit = 1000 !< limit of file (for one item)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: FILE_EXTERNAL_INPUT_time_advance

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: I_prev = 1 !< [index] previous
  integer, private, parameter :: I_next = 2 !< [index] next

  integer, private, parameter :: I_periodic_year  = 1
  integer, private, parameter :: I_periodic_month = 2
  integer, private, parameter :: I_periodic_day   = 3

  integer, private, parameter :: FILE_EXTERNAL_INPUT_item_limit = 1000  !< limit of item
  integer, private, parameter :: FILE_EXTERNAL_INPUT_step_limit = 10000 !< limit of steps          for each item
  integer, private, parameter :: FILE_EXTERNAL_INPUT_dim_limit  = 3     !< limit of dimension rank for each item

  type, private :: itemcontainer
     character(len=H_SHORT)             :: varname                   !< variable name
     integer                            :: nfile                     !< number of files
     integer                            :: file_current              !< current number of the file
     character(len=H_LONG), allocatable :: basename(:)               !< file names
     integer                            :: fid                       !< file id
     integer                            :: ndim                      !< number of dimensions
     integer                            :: dim_size(FILE_EXTERNAL_INPUT_dim_limit) !< size of dimension (z,x,y)
     integer, allocatable               :: dim_start(:)              !< start index
     integer                            :: step_limit                !< size limit of time dimension
     integer                            :: step_num                  !< size of time dimension
     real(DP), allocatable              :: time(:)                   !< time of each step [sec]
     logical                            :: fixed_step                !< fix step position?
     integer                            :: flag_periodic             !< treat as periodic data? (0:no 1:yearly 2:monthly 3:daily)
     real(RP)                           :: offset                    !< offset value (default is set to 0)
     integer                            :: data_step_prev            !< step position to read, previous from current
     integer                            :: data_step_next            !< step position to read, next     to   current
     integer                            :: data_step_offset          !< offset of step position for each file
     real(RP), allocatable              :: value(:,:,:,:)            !< data value                    (1:previous,2:next)
     character(len=H_SHORT)             :: axistype                  !< axis type
     logical                            :: transpose                 !< true: xyz, false: zxy
  end type itemcontainer

  integer,             private :: FILE_EXTERNAL_INPUT_item_count = 0                       !< number of item to output
  type(itemcontainer), private :: FILE_EXTERNAL_INPUT_item(FILE_EXTERNAL_INPUT_item_limit) !< item to output

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILE_EXTERNAL_INPUT_setup
    use scale_process, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=H_LONG)  :: basename(FILE_EXTERNAL_INPUT_file_limit)
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

    namelist /EXTERNAL_ITEM/ &
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
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[EXTERNAL_INPUT] / Categ[FILE] / Origin[SCALElib]'

    ! count external data from namelist
    rewind(IO_FID_CONF)
    do count = 1, FILE_EXTERNAL_INPUT_item_limit
       ! set default
       step_limit            = FILE_EXTERNAL_INPUT_step_limit
       basename(:)           = ''
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
       read(IO_FID_CONF,nml=EXTERNAL_ITEM,iostat=ierr)
       if ( ierr < 0 ) then !--- no more items
          exit
       elseif( ierr > 0 ) then !--- fatal error
          write(*,*) 'xxx Not appropriate names in namelist EXTERNAL_ITEM. Check!', count
          call PRC_abort
       endif
       if( IO_NML .AND. IO_FID_NML /= IO_FID_LOG ) write(IO_FID_NML,nml=EXTERNAL_ITEM)

       call FILE_EXTERNAL_INPUT_regist( basename(:),           & ! [IN]
                          varname,               & ! [IN]
                          axistype,              & ! [IN]
                          enable_periodic_year,  & ! [IN]
                          enable_periodic_month, & ! [IN]
                          enable_periodic_day,   & ! [IN]
                          step_fixed,            & ! [IN]
                          offset,                & ! [IN]
                          defval,                & ! [IN]
                          check_coordinates,     & ! [IN]
                          step_limit             ) ! [IN]
    enddo

    return
  end subroutine FILE_EXTERNAL_INPUT_setup

  !-----------------------------------------------------------------------------
  !> Regist data
  subroutine FILE_EXTERNAL_INPUT_regist( &
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
       step_limit,            &
       exist                  )
    use scale_file_h, only: &
       FILE_FREAD
    use scale_file, only: &
       FILE_Open,             &
       FILE_Get_All_DataInfo, &
       FILE_Read
    use scale_process, only: &
       PRC_myrank, &
       PRC_abort
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
       TIME_NOWDAYSEC,   &
       TIME_OFFSET_year
    use scale_file_cartesC, only: &
       FILE_CARTESC_check_coordinates
    implicit none

    character(len=*), intent(in)  :: basename(FILE_EXTERNAL_INPUT_file_limit)
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype
    integer,          intent(in)  :: step_fixed            ! fixed step position to read
    logical,          intent(in)  :: enable_periodic_year  ! treat as yearly               periodic data?
    logical,          intent(in)  :: enable_periodic_month ! treat as yearly,monthly       periodic data?
    logical,          intent(in)  :: enable_periodic_day   ! treat as yearly,monthly,daily periodic data?
    real(RP),         intent(in)  :: offset
    real(RP),         intent(in)  :: defval

    logical,          intent(in),  optional :: check_coordinates
    integer,          intent(in),  optional :: step_limit            ! limit number for reading data
    logical,          intent(out), optional :: exist

    integer                :: step_nmax
    character(len=H_LONG)  :: description
    character(len=H_SHORT) :: unit
    integer                :: datatype
    integer                :: dim_rank
    character(len=H_SHORT) :: dim_name  (3)
    integer                :: dim_size  (3)
    real(DP)               :: time_start(FILE_EXTERNAL_INPUT_step_limit)
    real(DP)               :: time_end  (FILE_EXTERNAL_INPUT_step_limit)
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
    !---------------------------------------------------------------------------

    if ( present(step_limit) ) then
       if ( step_limit > 0 ) then
          step_limit_ = step_limit
       else
          step_limit_ = FILE_EXTERNAL_INPUT_step_limit
       endif
    else
       step_limit_ = FILE_EXTERNAL_INPUT_step_limit
    endif

    do nid = 1, FILE_EXTERNAL_INPUT_item_count
       if ( FILE_EXTERNAL_INPUT_item(nid)%varname  == varname ) then
          write(*,*) 'xxx Data is already registered! basename,varname = ', trim(basename(1)), ', ', trim(varname)
          call PRC_abort
       endif
    enddo

    FILE_EXTERNAL_INPUT_item_count = FILE_EXTERNAL_INPUT_item_count + 1

    if ( FILE_EXTERNAL_INPUT_item_count > FILE_EXTERNAL_INPUT_item_limit ) then
       write(*,*) 'xxx Number of EXT data exceedes the limit', FILE_EXTERNAL_INPUT_item_count, FILE_EXTERNAL_INPUT_item_limit
       call PRC_abort
    endif

    call FILE_Open( basename(1),      & ! [IN]
                    fid,              & ! [OUT]
                    rankid=PRC_myrank ) ! [IN]

    ! read from file
    call FILE_Get_All_Datainfo( step_limit_,                   & ! [IN]
                                FILE_EXTERNAL_INPUT_dim_limit, & ! [IN]
                                fid,                           & ! [IN]
                                varname,                       & ! [IN]
                                step_nmax,                     & ! [OUT]
                                description,                   & ! [OUT]
                                unit,                          & ! [OUT]
                                datatype,                      & ! [OUT]
                                dim_rank,                      & ! [OUT]
                                dim_name  (:),                 & ! [OUT]
                                dim_size  (:),                 & ! [OUT]
                                time_start(1:step_limit_),     & ! [OUT]
                                time_end  (1:step_limit_),     & ! [OUT]
                                time_units                     ) ! [OUT]

    if ( step_nmax > 0 ) then
       if ( present(exist) ) then
          exist = .true.
       endif
    else
       if ( present(exist) ) then
          exist = .false.
          return
       else
          write(*,*) 'xxx Data not found! basename,varname = ', trim(basename(1)), ', ', trim(varname)
          call PRC_abort
       endif
    endif

    do n = dim_rank+1, 3
       dim_size(n) = 1
    enddo

    nid = FILE_EXTERNAL_INPUT_item_count

    do n = 1, FILE_EXTERNAL_INPUT_file_limit
       if( basename(n) == '' ) exit
    enddo
    FILE_EXTERNAL_INPUT_item(nid)%nfile            = n - 1
    FILE_EXTERNAL_INPUT_item(nid)%file_current     = 1
    FILE_EXTERNAL_INPUT_item(nid)%data_step_offset = 0

    allocate( FILE_EXTERNAL_INPUT_item(nid)%basename(FILE_EXTERNAL_INPUT_item(nid)%nfile) )
    FILE_EXTERNAL_INPUT_item(nid)%basename(1:FILE_EXTERNAL_INPUT_item(nid)%nfile) = basename(1:FILE_EXTERNAL_INPUT_item(nid)%nfile)

    ! setup item
    FILE_EXTERNAL_INPUT_item(nid)%fid         = fid
    FILE_EXTERNAL_INPUT_item(nid)%varname     = varname
    FILE_EXTERNAL_INPUT_item(nid)%dim_size(:) = dim_size(:)
    FILE_EXTERNAL_INPUT_item(nid)%step_num    = step_nmax
    FILE_EXTERNAL_INPUT_item(nid)%step_limit  = step_limit_

    if ( enable_periodic_day ) then
       FILE_EXTERNAL_INPUT_item(nid)%flag_periodic = I_periodic_day
    elseif( enable_periodic_month ) then
       FILE_EXTERNAL_INPUT_item(nid)%flag_periodic = I_periodic_month
    elseif( enable_periodic_year ) then
       FILE_EXTERNAL_INPUT_item(nid)%flag_periodic = I_periodic_year
    else
       FILE_EXTERNAL_INPUT_item(nid)%flag_periodic = 0
    endif

    allocate( FILE_EXTERNAL_INPUT_item(nid)%value(dim_size(1),dim_size(2),dim_size(3),2) )
    FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,:) = defval
    FILE_EXTERNAL_INPUT_item(nid)%offset         = offset

    allocate( FILE_EXTERNAL_INPUT_item(nid)%time(step_limit_) )
    FILE_EXTERNAL_INPUT_item(nid)%time(:) = 0.0_DP

    do n = 1, FILE_EXTERNAL_INPUT_item(nid)%step_num
       FILE_EXTERNAL_INPUT_item(nid)%time(n) = CALENDAR_CFunits2sec( time_end(n), time_units, TIME_OFFSET_year, TIME_STARTDAYSEC )
    enddo

    if ( FILE_EXTERNAL_INPUT_item(nid)%step_num == 1 ) then

       FILE_EXTERNAL_INPUT_item(nid)%fixed_step     = .true.
       FILE_EXTERNAL_INPUT_item(nid)%data_step_prev = 1
       FILE_EXTERNAL_INPUT_item(nid)%data_step_next = 1

    else if ( step_fixed > 0 ) then ! fixed time step mode

       FILE_EXTERNAL_INPUT_item(nid)%fixed_step     = .true.
       FILE_EXTERNAL_INPUT_item(nid)%data_step_prev = step_fixed
       FILE_EXTERNAL_INPUT_item(nid)%data_step_next = step_fixed

    else

       FILE_EXTERNAL_INPUT_item(nid)%fixed_step = .false.

       ! seek start position
       FILE_EXTERNAL_INPUT_item(nid)%data_step_next = 1
       do n = 1, FILE_EXTERNAL_INPUT_item(nid)%step_num
          if ( FILE_EXTERNAL_INPUT_item(nid)%time(n) > TIME_NOWDAYSEC ) exit
          FILE_EXTERNAL_INPUT_item(nid)%data_step_next = n + 1
       enddo

       FILE_EXTERNAL_INPUT_item(nid)%data_step_prev = FILE_EXTERNAL_INPUT_item(nid)%data_step_next - 1

       if ( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic > 0 ) then ! periodic time step mode

          if ( FILE_EXTERNAL_INPUT_item(nid)%data_step_next == 1 ) then ! between first-1 and first

             ! first-1 = last
             FILE_EXTERNAL_INPUT_item(nid)%data_step_prev = FILE_EXTERNAL_INPUT_item(nid)%step_num

          elseif( FILE_EXTERNAL_INPUT_item(nid)%data_step_next == FILE_EXTERNAL_INPUT_item(nid)%step_num+1 ) then ! between last and last+1

             ! last+1 = first
             FILE_EXTERNAL_INPUT_item(nid)%data_step_next = 1

             ! update data time in periodic condition
             do n = 1, FILE_EXTERNAL_INPUT_item(nid)%step_num
                dataday     = 0
                datasec     = FILE_EXTERNAL_INPUT_item(nid)%time(n)
                offset_year = 0
                call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

                call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                           datasubsec,  & ! [OUT]
                                           dataday,     & ! [IN]
                                           datasec,     & ! [IN]
                                           offset_year  ) ! [IN]

                if    ( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_day   ) then
                   datadate(I_day)   = datadate(I_day)   + 1
                elseif( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_month ) then
                   datadate(I_month) = datadate(I_month) + 1
                elseif( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_year  ) then
                   datadate(I_year)  = datadate(I_year)  + 1
                endif

                call CALENDAR_date2daysec( dataday,     & ! [OUT]
                                           datasec,     & ! [OUT]
                                           datadate(:), & ! [IN]
                                           datasubsec,  & ! [IN]
                                           offset_year  ) ! [IN]

                FILE_EXTERNAL_INPUT_item(nid)%time(n) = CALENDAR_combine_daysec( dataday, datasec )
             enddo

             if( IO_L ) write(IO_FID_LOG,*) '*** data time is updated.'
          endif

       else ! normal mode

          if (      FILE_EXTERNAL_INPUT_item(nid)%data_step_next == 1                          &
               .OR. FILE_EXTERNAL_INPUT_item(nid)%data_step_next == FILE_EXTERNAL_INPUT_item(nid)%step_num+1 ) then
             write(*,*) 'xxx Current time is out of period of external data! ', trim(varname)
             call PRC_abort
          endif

       endif

    endif

    !--- read first data
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Initial read of external data : ', trim(varname)

    if (       dim_size(1) >= 1 &
         .AND. dim_size(2) == 1 &
         .AND. dim_size(3) == 1 ) then ! 1D

       call FILE_EXTERNAL_INPUT_get_dims1D( dim1_max, dim1_S, dim1_E, & ! [OUT]
                                            varname, axistype         ) ! [IN]

       FILE_EXTERNAL_INPUT_item(nid)%ndim         = 1
       FILE_EXTERNAL_INPUT_item(nid)%transpose    = .false.
       allocate( FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) )
       FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) = dim1_S

       if ( dim1_max /= dim_size(1) ) then
          write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
          write(*,*) 'xxx dim 1 (data,requested)    : ', dim_size(1), dim1_max
          call PRC_abort
       endif

       ! read prev
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A)') &
                  '*** Read 1D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_prev, ')'

       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,1,1,I_prev), & ! [OUT]
                       step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev  ) ! [IN]
       ! read next
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A)') &
                  '*** Read 1D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ')'

       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,1,1,I_next), & ! [OUT]
                       step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next  ) ! [IN]

    elseif(       dim_size(1) >= 1 &
            .AND. dim_size(2) >  1 &
            .AND. dim_size(3) == 1 ) then ! 2D

       call FILE_EXTERNAL_INPUT_get_dims2D( dim1_max, dim1_S, dim1_E,                & ! [OUT]
                                            dim2_max, dim2_S, dim2_E,                & ! [OUT]
                                            FILE_EXTERNAL_INPUT_item(nid)%transpose, & ! [OUT]
                                            varname, axistype                        ) ! [IN]

       FILE_EXTERNAL_INPUT_item(nid)%ndim         = 2
       allocate( FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) )
       FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) = dim1_S
       FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) = dim2_S

       if (      dim1_max /= dim_size(1) &
            .OR. dim2_max /= dim_size(2) ) then
          write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
          write(*,*) 'xxx dim 1 (data,requested)    : ', dim_size(1), dim1_max
          write(*,*) 'xxx dim 2 (data,requested)    : ', dim_size(2), dim2_max
          call PRC_abort
       endif

       ! read prev
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A)') &
                  '*** Read 2D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_prev, ')'

       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,1,I_prev), & ! [OUT]
                       step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev  ) ! [IN]
       ! read next
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A)') &
                  '*** Read 2D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ')'

       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,1,I_next), & ! [OUT]
                       step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next  ) ! [IN]

    elseif(       dim_size(1) >= 1 &
            .AND. dim_size(2) >  1 &
            .AND. dim_size(3) >  1 ) then ! 3D

       call FILE_EXTERNAL_INPUT_get_dims3D( dim1_max, dim1_S, dim1_E,                & ! [OUT]
                                            dim2_max, dim2_S, dim2_E,                & ! [OUT]
                                            dim3_max, dim3_S, dim3_E,                & ! [OUT]
                                            FILE_EXTERNAL_INPUT_item(nid)%transpose, & ! [OUT]
                                            varname, axistype                        ) ! [IN]

       FILE_EXTERNAL_INPUT_item(nid)%ndim         = 3
       allocate( FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) )
       FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) = dim1_S
       FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) = dim2_S
       FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) = dim3_S

       if (      dim1_max /= dim_size(1) &
            .OR. dim2_max /= dim_size(2) &
            .OR. dim3_max /= dim_size(3) ) then
          write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
          write(*,*) 'xxx dim 1 (data,requested)    : ', dim_size(1), dim1_max
          write(*,*) 'xxx dim 2 (data,requested)    : ', dim_size(2), dim2_max
          write(*,*) 'xxx dim 3 (data,requested)    : ', dim_size(3), dim3_max
          call PRC_abort
       endif

       ! read prev
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A)') &
                  '*** Read 3D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_prev, ')'

       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev), & ! [OUT]
                       step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev  ) ! [IN]

       ! read next
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A)') &
                  '*** Read 3D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ')'

       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next), & ! [OUT]
                       step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next  ) ! [IN]

    else
       write(*,*) 'xxx Unexpected dimsize: ', dim_size(:)
       call PRC_abort
    endif

    if ( present(check_coordinates) ) then
       if ( check_coordinates ) then
          call FILE_CARTESC_check_coordinates( fid,                                                &
                                               atmos     = FILE_EXTERNAL_INPUT_item(nid)%ndim==3,  &
                                               transpose = FILE_EXTERNAL_INPUT_item(nid)%transpose )
       endif
    endif

    return
  end subroutine FILE_EXTERNAL_INPUT_regist

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine FILE_EXTERNAL_INPUT_update_1D( &
       varname,      &
       time_current, &
       var,          &
       error         )
    use scale_file, only: &
       FILE_Read
    use scale_process, only: &
       PRC_abort
    implicit none
    character(len=*), intent(in)  :: varname      ! item name
    real(DP),         intent(in)  :: time_current ! current time
    real(RP),         intent(out) :: var(:)       ! variable
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile
    integer  :: step_next

    integer  :: n
    integer  :: n1
    integer  :: nn1
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, FILE_EXTERNAL_INPUT_item_count
       if( varname == FILE_EXTERNAL_INPUT_item(n)%varname ) nid = n
    enddo

    if ( nid == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Variable was not registered: ', trim(varname)
       return
    endif

    if ( FILE_EXTERNAL_INPUT_item(nid)%ndim /= 1 ) then
       write(*,*) 'xxx Data is not 1D var: ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
       call PRC_abort
    endif

    call FILE_EXTERNAL_INPUT_time_advance( nid,          & ! [IN]
                                           time_current, & ! [IN]
                                           weight,       & ! [OUT]
                                           do_readfile   ) ! [OUT]

    if ( do_readfile ) then
       step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next - FILE_EXTERNAL_INPUT_item(nid)%data_step_offset

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A,I4,A)') &
                  '*** Read 1D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ', file step=', step_next, ')'

       ! next -> prev
       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

       ! read next
       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,1,1,I_next), & ! [OUT]
                       step=step_next                       ) ! [IN]
    endif

    ! store data with weight
    do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(1)
       nn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1

       var(nn1) = ( 1.0_RP-weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,1,1,I_prev) &
                + (        weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,1,1,I_next)
    enddo

    error = .false.

    return
  end subroutine FILE_EXTERNAL_INPUT_update_1D

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine FILE_EXTERNAL_INPUT_update_2D( &
       varname,      &
       time_current, &
       var,          &
       error         )
    use scale_file, only: &
       FILE_Read
    implicit none
    character(len=*), intent(in)  :: varname      ! item name
    real(DP),         intent(in)  :: time_current ! current time
    real(RP),         intent(out) :: var(:,:)     ! variable
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile
    integer  :: step_next

    integer  :: n
    integer  :: n1, n2
    integer  :: nn1, nn2
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, FILE_EXTERNAL_INPUT_item_count
       if( varname == FILE_EXTERNAL_INPUT_item(n)%varname ) nid = n
    enddo

    if ( nid == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx variable was not registered: ', trim(varname)
       return
    endif

    call FILE_EXTERNAL_INPUT_time_advance( nid,          & ! [IN]
                                           time_current, & ! [IN]
                                           weight,       & ! [OUT]
                                           do_readfile   ) ! [OUT]

    if ( do_readfile ) then

       step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next - FILE_EXTERNAL_INPUT_item(nid)%data_step_offset

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A,I4,A)') &
                  '*** Read 2D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ', file step=', step_next, ')'

       ! next -> prev
       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

       ! read next
       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,1,I_next), & ! [OUT]
                       step=step_next                       ) ! [IN]
    endif

    if ( FILE_EXTERNAL_INPUT_item(nid)%transpose ) then
       ! store data with weight (x,z)->(z,x)
       do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(1)
          nn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1

          do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(2)
          nn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1

             var(nn2,nn1) = ( 1.0_RP-weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,n2,1,I_prev) &
                          + (        weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,n2,1,I_next)
          enddo
       enddo
    else
       ! store data with weight
       do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(2)
          nn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1

          do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(1)
             nn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1

             var(nn1,nn2) = ( 1.0_RP-weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,n2,1,I_prev) &
                          + (        weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,n2,1,I_next)
          enddo
       enddo
    endif

    error = .false.

    return
  end subroutine FILE_EXTERNAL_INPUT_update_2D

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine FILE_EXTERNAL_INPUT_update_3D( &
       varname,      &
       time_current, &
       var,          &
       error         )
    use scale_file, only: &
       FILE_Read
    implicit none
    character(len=*), intent(in)  :: varname      ! item name
    real(DP),         intent(in)  :: time_current ! current time
    real(RP),         intent(out) :: var(:,:,:)   ! variable
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile
    integer  :: step_next

    integer  :: n
    integer  :: n1, n2, n3
    integer  :: nn1, nn2, nn3
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, FILE_EXTERNAL_INPUT_item_count
       if( varname == FILE_EXTERNAL_INPUT_item(n)%varname ) nid = n
    enddo

    if ( nid == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx variable was not registered: ', trim(varname)
       return
    endif

    call FILE_EXTERNAL_INPUT_time_advance( nid,          & ! [IN]
                                           time_current, & ! [IN]
                                           weight,       & ! [OUT]
                                           do_readfile   ) ! [OUT]

    if ( do_readfile ) then

       step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next - FILE_EXTERNAL_INPUT_item(nid)%data_step_offset

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I4,A,I4,A)') &
                  '*** Read 3D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ', file step=', step_next, ')'

       ! next -> prev
       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

       ! read next
       call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                 & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                       FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next), & ! [OUT]
                       step=step_next                       ) ! [IN]
    endif

    if ( FILE_EXTERNAL_INPUT_item(nid)%transpose ) then
       ! store data with weight (x,y,z)->(z,x,y)
       do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(2)
          nn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1

          do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(1)
             nn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1

             do n3 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(3)
                nn3 = n3 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) - 1

                var(nn3,nn1,nn2) = ( 1.0_RP-weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,n2,n3,I_prev) &
                                 + (        weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,n2,n3,I_next)
             enddo
          enddo
       enddo
    else
       ! store data with weight (z,x,y)->(z,x,y)
       do n3 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(3)
          nn3 = n3 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) - 1

          do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(2)
             nn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1

             do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_size(1)
                nn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1

                var(nn1,nn2,nn3) = ( 1.0_RP-weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,n2,n3,I_prev) &
                                 + (        weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(n1,n2,n3,I_next)
             enddo
          enddo
       enddo
    endif

    error = .false.

    return
  end subroutine FILE_EXTERNAL_INPUT_update_3D

  !-----------------------------------------------------------------------------
  !> Check time to read and calc weight
  subroutine FILE_EXTERNAL_INPUT_time_advance( &
       nid,          &
       time_current, &
       weight,       &
       do_readfile   )
    use scale_file_h, only: &
       FILE_FREAD
    use scale_file, only: &
       FILE_Open,           &
       FILE_Get_All_DataInfo
    use scale_process, only: &
       PRC_myrank, &
       PRC_abort
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
       TIME_OFFSET_year
    implicit none

    integer,  intent(in)  :: nid          ! item id
    real(DP), intent(in)  :: time_current ! current time
    real(RP), intent(out) :: weight       ! weight
    logical,  intent(out) :: do_readfile  ! read new data at this time?

    integer                :: step_nmax
    character(len=H_LONG)  :: description
    character(len=H_SHORT) :: unit
    integer                :: datatype
    integer                :: dim_rank
    character(len=H_SHORT) :: dim_name  (3)
    integer                :: dim_size  (3)
    real(DP)               :: time_start(FILE_EXTERNAL_INPUT_step_limit)
    real(DP)               :: time_end  (FILE_EXTERNAL_INPUT_step_limit)
    character(len=H_MID)   :: time_units

    integer  :: datadate(6)   !< date
    real(DP) :: datasubsec    !< subsecond
    integer  :: dataday       !< absolute day
    real(DP) :: datasec       !< absolute second
    integer  :: offset_year   !< offset year

    real(DP) :: time_prev, time_next
    integer  :: step_prev, step_next
    integer  :: t
    integer  :: fid
    integer  :: n, nn
    !---------------------------------------------------------------------------

    do_readfile = .false.

    if ( FILE_EXTERNAL_INPUT_item(nid)%fixed_step ) then
       !--- no time-advance
    else
       ! time is passed?
       if ( time_current > FILE_EXTERNAL_INPUT_item(nid)%time( FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ) then

          do_readfile = .true.

          if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Update external input : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)

          ! update step position
          FILE_EXTERNAL_INPUT_item(nid)%data_step_prev = FILE_EXTERNAL_INPUT_item(nid)%data_step_prev + 1
          FILE_EXTERNAL_INPUT_item(nid)%data_step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next + 1

          if ( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic > 0 ) then ! periodic time step mode

             if ( FILE_EXTERNAL_INPUT_item(nid)%data_step_next == FILE_EXTERNAL_INPUT_item(nid)%step_num+1 ) then

                ! last+1 = first
                FILE_EXTERNAL_INPUT_item(nid)%data_step_next = 1

                ! update data time in periodic condition
                do t = 1, FILE_EXTERNAL_INPUT_item(nid)%step_num
                   dataday     = 0
                   datasec     = FILE_EXTERNAL_INPUT_item(nid)%time(t)
                   offset_year = 0
                   call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

                   call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                              datasubsec,  & ! [OUT]
                                              dataday,     & ! [IN]
                                              datasec,     & ! [IN]
                                              offset_year  ) ! [IN]

                   if    ( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_day   ) then
                      datadate(I_day)   = datadate(I_day)   + 1
                   elseif( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_month ) then
                      datadate(I_month) = datadate(I_month) + 1
                   elseif( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_year  ) then
                      datadate(I_year)  = datadate(I_year)  + 1
                   endif

                   call CALENDAR_date2daysec( dataday,     & ! [IN]
                                              datasec,     & ! [IN]
                                              datadate(:), & ! [OUT]
                                              datasubsec,  & ! [OUT]
                                              offset_year  ) ! [IN]

                   FILE_EXTERNAL_INPUT_item(nid)%time(t) = CALENDAR_combine_daysec( dataday, datasec )
                enddo
             endif

          else ! normal mode

             if ( FILE_EXTERNAL_INPUT_item(nid)%data_step_next == FILE_EXTERNAL_INPUT_item(nid)%step_num+1 ) then

                if ( FILE_EXTERNAL_INPUT_item(nid)%file_current < FILE_EXTERNAL_INPUT_item(nid)%nfile ) then

                   FILE_EXTERNAL_INPUT_item(nid)%file_current = FILE_EXTERNAL_INPUT_item(nid)%file_current + 1

                   call FILE_Open( FILE_EXTERNAL_INPUT_item(nid)%basename(FILE_EXTERNAL_INPUT_item(nid)%file_current), & ! [IN]
                                   fid,              & ! [OUT]
                                   rankid=PRC_myrank ) ! [IN]

                   ! read from file
                   call FILE_Get_All_Datainfo( FILE_EXTERNAL_INPUT_item(nid)%step_limit,               & ! [IN]
                                               FILE_EXTERNAL_INPUT_dim_limit,                          & ! [IN]
                                               fid,                                                    & ! [IN]
                                               FILE_EXTERNAL_INPUT_item(nid)%varname,                  & ! [IN]
                                               step_nmax,                                              & ! [OUT]
                                               description,                                            & ! [OUT]
                                               unit,                                                   & ! [OUT]
                                               datatype,                                               & ! [OUT]
                                               dim_rank,                                               & ! [OUT]
                                               dim_name  (:),                                          & ! [OUT]
                                               dim_size  (:),                                          & ! [OUT]
                                               time_start(1:FILE_EXTERNAL_INPUT_item(nid)%step_limit), & ! [OUT]
                                               time_end  (1:FILE_EXTERNAL_INPUT_item(nid)%step_limit), & ! [OUT]
                                               time_units                                              ) ! [OUT]

                   if ( step_nmax == 0 ) then
                      write(*,*) 'xxx Data not found! basename = ', trim(FILE_EXTERNAL_INPUT_item(nid)%basename(FILE_EXTERNAL_INPUT_item(nid)%file_current)), &
                                                    ', varname = ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
                      call PRC_abort
                   endif

                   do n = 1, dim_rank
                      if ( FILE_EXTERNAL_INPUT_item(nid)%dim_size(n) /= dim_size(n) ) then
                         write(*,*) 'xxx The size of dimension', n, ' is inconsistent! '
                         write(*,*) 'xxx size (previous,current) = ', FILE_EXTERNAL_INPUT_item(nid)%dim_size(n), dim_size(n)
                         write(*,*) 'xxx basename = ', trim(FILE_EXTERNAL_INPUT_item(nid)%basename(FILE_EXTERNAL_INPUT_item(nid)%file_current)), &
                                       ', varname = ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
                         call PRC_abort
                      endif
                   enddo

                   do n = 1, step_nmax
                      nn = FILE_EXTERNAL_INPUT_item(nid)%step_num + n
                      FILE_EXTERNAL_INPUT_item(nid)%time(nn) = CALENDAR_CFunits2sec( time_end(n), time_units, TIME_OFFSET_year, TIME_STARTDAYSEC )
                   enddo

                   if ( FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_prev) > FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_next) ) then
                      write(*,*) 'xxx Time in new file is earlier than last time of previous file! stop'
                      write(*,*) 'xxx Time (previous,current)  = ', FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_prev), &
                                                                    FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_next)
                      write(*,*) 'xxx Data not found! basename = ', trim(FILE_EXTERNAL_INPUT_item(nid)%basename(FILE_EXTERNAL_INPUT_item(nid)%file_current)), &
                                                    ', varname = ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
                      call PRC_abort
                   endif

                   FILE_EXTERNAL_INPUT_item(nid)%fid              = fid
                   FILE_EXTERNAL_INPUT_item(nid)%data_step_offset = FILE_EXTERNAL_INPUT_item(nid)%step_num
                   FILE_EXTERNAL_INPUT_item(nid)%step_num         = FILE_EXTERNAL_INPUT_item(nid)%step_num + step_nmax

                else
                   write(*,*) 'xxx Current time is out of period of external data! '
                   call PRC_abort
                endif

             endif

          endif ! periodic or not

       endif ! time is passed?

    endif ! fixed step or not

    ! calc weight
    if ( FILE_EXTERNAL_INPUT_item(nid)%fixed_step ) then

       weight = 0.0_RP

    elseif( FILE_EXTERNAL_INPUT_item(nid)%data_step_next == 1 ) then ! periodic case

       step_prev = FILE_EXTERNAL_INPUT_item(nid)%data_step_prev
       step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next

       dataday     = 0
       datasec     = FILE_EXTERNAL_INPUT_item(nid)%time( step_prev )
       offset_year = 0
       call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

       call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                  datasubsec,  & ! [OUT]
                                  dataday,     & ! [IN]
                                  datasec,     & ! [IN]
                                  offset_year  ) ! [IN]

       if    ( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_day   ) then
          datadate(I_day)   = datadate(I_day)   - 1
       elseif( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_month ) then
          datadate(I_month) = datadate(I_month) - 1
       elseif( FILE_EXTERNAL_INPUT_item(nid)%flag_periodic == I_periodic_year  ) then
          datadate(I_year)  = datadate(I_year)  - 1
       endif

       call CALENDAR_date2daysec( dataday,     & ! [IN]
                                  datasec,     & ! [IN]
                                  datadate(:), & ! [OUT]
                                  datasubsec,  & ! [OUT]
                                  offset_year  ) ! [IN]

       time_prev = CALENDAR_combine_daysec( dataday, datasec )
       time_next = FILE_EXTERNAL_INPUT_item(nid)%time( step_next )

       weight = ( time_current - time_prev ) &
              / ( time_next    - time_prev )

    else ! normal case

       step_prev = FILE_EXTERNAL_INPUT_item(nid)%data_step_prev
       step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next

       time_prev = FILE_EXTERNAL_INPUT_item(nid)%time( step_prev )
       time_next = FILE_EXTERNAL_INPUT_item(nid)%time( step_next )

       weight = ( time_current - time_prev ) &
              / ( time_next    - time_prev )

    endif

    return
  end subroutine FILE_EXTERNAL_INPUT_time_advance

end module scale_file_external_input
