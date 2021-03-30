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
#include "scalelib.h"
module scale_file_external_input
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
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
  public :: FILE_EXTERNAL_INPUT_put_ref
  public :: FILE_EXTERNAL_INPUT_get_ref
  public :: FILE_EXTERNAL_INPUT_query

  interface FILE_EXTERNAL_INPUT_regist
     module procedure FILE_EXTERNAL_INPUT_regist_file
     module procedure FILE_EXTERNAL_INPUT_regist_external_1D
     module procedure FILE_EXTERNAL_INPUT_regist_external_2D
     module procedure FILE_EXTERNAL_INPUT_regist_external_3D
  end interface FILE_EXTERNAL_INPUT_regist

  interface FILE_EXTERNAL_INPUT_update
     module procedure FILE_EXTERNAL_INPUT_update_1D
     module procedure FILE_EXTERNAL_INPUT_update_2D
     module procedure FILE_EXTERNAL_INPUT_update_3D
  end interface FILE_EXTERNAL_INPUT_update

  interface FILE_EXTERNAL_INPUT_put_ref
     module procedure FILE_EXTERNAL_INPUT_put_ref_1D
     module procedure FILE_EXTERNAL_INPUT_put_ref_2D
     module procedure FILE_EXTERNAL_INPUT_put_ref_3D
  end interface FILE_EXTERNAL_INPUT_put_ref

  interface FILE_EXTERNAL_INPUT_get_ref
     module procedure FILE_EXTERNAL_INPUT_get_ref_1D
     module procedure FILE_EXTERNAL_INPUT_get_ref_2D
     module procedure FILE_EXTERNAL_INPUT_get_ref_3D
  end interface FILE_EXTERNAL_INPUT_get_ref

  abstract interface
     subroutine get_dims1D( &
          dim1_size, &
          dim1_max,  &
          dim1_S,    &
          varname,   &
          axistype   )
       integer,          intent(out) :: dim1_size
       integer,          intent(out) :: dim1_max
       integer,          intent(out) :: dim1_S
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: axistype     ! axis type (Z/X/Y)
     end subroutine get_dims1D

     subroutine get_dims2D( &
          dim1_size, &
          dim1_max,  &
          dim1_S,    &
          dim2_size, &
          dim2_max,  &
          dim2_S,    &
          transpose, &
          varname,   &
          axistype   )
       integer,          intent(out) :: dim1_size
       integer,          intent(out) :: dim1_max
       integer,          intent(out) :: dim1_S
       integer,          intent(out) :: dim2_size
       integer,          intent(out) :: dim2_max
       integer,          intent(out) :: dim2_S
       logical,          intent(out) :: transpose
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: axistype     ! axis type (XY/XZ/ZX)
     end subroutine get_dims2D

     subroutine get_dims3D( &
          dim1_size, &
          dim1_max,  &
          dim1_S,    &
          dim2_size, &
          dim2_max,  &
          dim2_S,    &
          dim3_size, &
          dim3_max,  &
          dim3_S,    &
          transpose, &
          varname,   &
          axistype   )
       integer,          intent(out) :: dim1_size
       integer,          intent(out) :: dim1_max
       integer,          intent(out) :: dim1_S
       integer,          intent(out) :: dim2_size
       integer,          intent(out) :: dim2_max
       integer,          intent(out) :: dim2_S
       integer,          intent(out) :: dim3_size
       integer,          intent(out) :: dim3_max
       integer,          intent(out) :: dim3_S
       logical,          intent(out) :: transpose
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: axistype     ! axis type (ZXY/XYZ/Land/Urban)
     end subroutine get_dims3D

     subroutine read1D( fid, varname, dim_type, &
                        var, &
                        step )
       use scale_precision
       integer,          intent(in)  :: fid
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: dim_type
       real(RP),         intent(out) :: var(:)
       integer,          intent(in), optional :: step
     end subroutine read1D

     subroutine read2D( fid, varname, dim_type, &
                         var, &
                         step )
       use scale_precision
       integer,          intent(in)  :: fid
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: dim_type
       real(RP),         intent(out) :: var(:,:)
       integer,          intent(in), optional :: step
     end subroutine read2D

     subroutine read3D( fid, varname, dim_type, &
                         var, &
                         step )
       use scale_precision
       integer,          intent(in)  :: fid
       character(len=*), intent(in)  :: varname
       character(len=*), intent(in)  :: dim_type
       real(RP),         intent(out) :: var(:,:,:)
       integer,          intent(in), optional :: step
     end subroutine read3D

  end interface

  procedure(get_dims1D), pointer :: FILE_EXTERNAL_INPUT_get_dims1D => NULL()
  procedure(get_dims2D), pointer :: FILE_EXTERNAL_INPUT_get_dims2D => NULL()
  procedure(get_dims3D), pointer :: FILE_EXTERNAL_INPUT_get_dims3D => NULL()
  public :: FILE_EXTERNAL_INPUT_get_dims1D
  public :: FILE_EXTERNAL_INPUT_get_dims2D
  public :: FILE_EXTERNAL_INPUT_get_dims3D

  procedure(read1D), pointer :: FILE_EXTERNAL_INPUT_read_1D => NULL()
  procedure(read2D), pointer :: FILE_EXTERNAL_INPUT_read_2D => NULL()
  procedure(read3D), pointer :: FILE_EXTERNAL_INPUT_read_3D => NULL()
  public :: FILE_EXTERNAL_INPUT_read_1D
  public :: FILE_EXTERNAL_INPUT_read_2D
  public :: FILE_EXTERNAL_INPUT_read_3D
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: I_prev = 1 !< [index] previous
  integer, public, parameter :: I_next = 2 !< [index] next

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: FILE_EXTERNAL_INPUT_query_Id
  private :: FILE_EXTERNAL_INPUT_time_advance
  private :: FILE_EXTERNAL_INPUT_init_var
  private :: FILE_EXTERNAL_INPUT_regist_var

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: I_periodic_year  = 1
  integer, private, parameter :: I_periodic_month = 2
  integer, private, parameter :: I_periodic_day   = 3

  integer, private, parameter :: FILE_EXTERNAL_INPUT_item_limit = 1000  !< limit of item
  integer, private, parameter :: FILE_EXTERNAL_INPUT_step_limit = 10000 !< limit of steps          for each item
  integer, private, parameter :: FILE_EXTERNAL_INPUT_dim_limit  = 3     !< limit of dimension rank for each item
  integer, private, parameter :: FILE_EXTERNAL_INPUT_att_limit  = 10    !< limit of dimension rank for each item

  type, private :: itemcontainer
     character(len=H_SHORT)             :: varname                   !< variable name
     logical                            :: file                      !< read from file?
     integer                            :: nfile                     !< number of files
     integer                            :: file_current              !< current number of the file
     character(len=H_LONG), allocatable :: basename(:)               !< file names
     integer                            :: fid                       !< file id
     integer                            :: ndim                      !< number of dimensions
     integer                            :: dim_size(FILE_EXTERNAL_INPUT_dim_limit) !< size of dimension (z,x,y)
     integer                            :: dim_start(FILE_EXTERNAL_INPUT_dim_limit)!< start index
     integer                            :: dim_max(FILE_EXTERNAL_INPUT_dim_limit)  !< number of valid data
     integer                            :: var_size(FILE_EXTERNAL_INPUT_dim_limit) !< size of input data
     integer                            :: var_start(FILE_EXTERNAL_INPUT_dim_limit)!< start index
     integer                            :: var_max(FILE_EXTERNAL_INPUT_dim_limit)  !< number of valid data
     integer                            :: step_limit                !< size limit of time dimension
     integer                            :: step_num                  !< size of time dimension
     real(DP), allocatable              :: time(:)                   !< time of each step [sec]
     logical                            :: fixed_step                !< fix step position?
     integer                            :: flag_periodic             !< treat as periodic data? (0:no 1:yearly 2:monthly 3:daily)
     integer                            :: data_step_prev            !< step position to read, previous from current
     integer                            :: data_step_next            !< step position to read, next     to   current
     integer                            :: data_step_offset          !< offset of step position for each file
     real(RP), allocatable              :: value(:,:,:,:)            !< data value                    (1:previous,2:next)
     character(len=H_SHORT)             :: axistype                  !< axis type
     logical                            :: transpose                 !< true: xyz, false: zxy
     logical                            :: aggregate                 !< file_aggregate
     logical                            :: allow_missing             !< error raised if false
  end type itemcontainer

  integer,             private :: FILE_EXTERNAL_INPUT_item_count = 0                       !< number of item to output
  type(itemcontainer), private :: FILE_EXTERNAL_INPUT_item(FILE_EXTERNAL_INPUT_item_limit) !< item to output

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILE_EXTERNAL_INPUT_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_file, only: &
       FILE_AGGREGATE_DEFAULT => FILE_AGGREGATE
    implicit none

    character(len=H_LONG)  :: basename
    logical                :: basename_add_num
    integer                :: number_of_files
    character(len=H_SHORT) :: varname
    character(len=H_SHORT) :: axistype
    integer                :: step_limit            ! limit number for reading data
    integer                :: step_fixed            ! fixed step position to read
    logical                :: enable_periodic_year  ! treat as yearly               periodic data?
    logical                :: enable_periodic_month ! treat as yearly,monthly       periodic data?
    logical                :: enable_periodic_day   ! treat as yearly,monthly,daily periodic data?
    real(RP)               :: defval
    logical                :: check_coordinates
    logical                :: file_aggregate
    logical                :: allow_missing

    namelist / EXTERNAL_ITEM / &
       basename,              &
       basename_add_num,      &
       number_of_files,       &
       varname,               &
       axistype,              &
       step_limit,            &
       step_fixed,            &
       enable_periodic_year,  &
       enable_periodic_month, &
       enable_periodic_day,   &
       defval,                &
       check_coordinates,     &
       file_aggregate,        &
       allow_missing

    integer  :: count
    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("FILE_EXTERNAL_INPUT_setup",*) 'Setup'

    ! count external data from namelist
    rewind(IO_FID_CONF)
    do count = 1, FILE_EXTERNAL_INPUT_item_limit
       ! set default
       step_limit            = FILE_EXTERNAL_INPUT_step_limit
       basename              = ''
       basename_add_num      = .false.
       number_of_files       = 1
       varname               = ''
       axistype              = ''
       step_fixed            = -1
       enable_periodic_year  = .false.
       enable_periodic_month = .false.
       enable_periodic_day   = .false.
       defval                = UNDEF
       check_coordinates     = .false.
       file_aggregate        = FILE_AGGREGATE_default
       allow_missing         = .false.

       ! read namelist
       read(IO_FID_CONF,nml=EXTERNAL_ITEM,iostat=ierr)
       if ( ierr < 0 ) then !--- no more items
          exit
       elseif( ierr > 0 ) then !--- fatal error
          LOG_ERROR("FILE_EXTERNAL_INPUT_setup",*) 'Not appropriate names in namelist EXTERNAL_ITEM. Check!', count
          call PRC_abort
       endif
       LOG_NML(EXTERNAL_ITEM)

       call FILE_EXTERNAL_INPUT_regist( basename,                              & ! [IN]
                                        basename_add_num,                      & ! [IN]
                                        number_of_files,                       & ! [IN]
                                        varname,                               & ! [IN]
                                        axistype,                              & ! [IN]
                                        enable_periodic_year,                  & ! [IN]
                                        enable_periodic_month,                 & ! [IN]
                                        enable_periodic_day,                   & ! [IN]
                                        step_fixed,                            & ! [IN]
                                        defval,                                & ! [IN]
                                        check_coordinates = check_coordinates, & ! [IN]
                                        aggregate         = file_aggregate,    & ! [IN]
                                        allow_missing     = allow_missing,     & ! [IN]
                                        step_limit        = step_limit         ) ! [IN]
    enddo

    return
  end subroutine FILE_EXTERNAL_INPUT_setup


  !-----------------------------------------------------------------------------
  !> check variable
  subroutine FILE_EXTERNAL_INPUT_init_var( &
       varname,       &
       axistype,      &
       aggregate,     &
       allow_missing  )
    use scale_prc, only: &
       PRC_abort
    implicit none
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype
    logical,          intent(in),  optional :: aggregate
    logical,          intent(in),  optional :: allow_missing

    logical :: aggregate_
    logical :: allow_missing_
    integer :: nid

    if ( present(aggregate) ) then
       aggregate_ = aggregate
    else
       aggregate_ = .false.
    end if

    if ( present(allow_missing) ) then
       allow_missing_ = allow_missing
    else
       allow_missing_ = .false.
    end if

    nid = FILE_EXTERNAL_INPUT_getId( varname )
    if ( nid > 0 ) then
       LOG_ERROR("FILE_EXTERNAL_INPUT_init_var",*) 'Data is already registered! varname = ', trim(varname)
       call PRC_abort
    endif

    FILE_EXTERNAL_INPUT_item_count = FILE_EXTERNAL_INPUT_item_count + 1

    if ( FILE_EXTERNAL_INPUT_item_count > FILE_EXTERNAL_INPUT_item_limit ) then
       LOG_ERROR("FILE_EXTERNAL_INPUT_init_var",*) 'Number of EXT data exceedes the limit', FILE_EXTERNAL_INPUT_item_count, FILE_EXTERNAL_INPUT_item_limit
       call PRC_abort
    endif

    nid = FILE_EXTERNAL_INPUT_item_count

    ! setup item
    FILE_EXTERNAL_INPUT_item(nid)%varname        = varname
    FILE_EXTERNAL_INPUT_item(nid)%axistype       = axistype
    FILE_EXTERNAL_INPUT_item(nid)%allow_missing  = allow_missing_
    FILE_EXTERNAL_INPUT_item(nid)%aggregate      = aggregate_
    FILE_EXTERNAL_INPUT_item(nid)%fixed_step     = .false.
    FILE_EXTERNAL_INPUT_item(nid)%flag_periodic  = 0
    FILE_EXTERNAL_INPUT_item(nid)%file           = .false.

    return
  end subroutine FILE_EXTERNAL_INPUT_init_var

  !-----------------------------------------------------------------------------
  !> Regist variable
  subroutine FILE_EXTERNAL_INPUT_regist_var( &
       varname,   &
       axistype,  &
       dim_rank,  &
       dim_size,  &
       step_num,  &
       file_num,  &
       time_now,  &
       time_step, &
       defval     )
    use scale_prc, only: &
       PRC_abort
    implicit none
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: axistype
    integer,          intent(in) :: dim_rank
    integer,          intent(in) :: dim_size(dim_rank)
    integer,          intent(in) :: step_num ! number of steps in single file
    integer,          intent(in) :: file_num ! number of files
    real(DP),         intent(in) :: time_now
    real(DP),         intent(in) :: time_step
    real(RP),         intent(in) :: defval

    integer  :: dim1_size, dim1_max, dim1_S
    integer  :: dim2_size, dim2_max, dim2_S
    integer  :: dim3_size, dim3_max, dim3_S

    integer :: nid
    integer :: n

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    select case ( dim_rank )
    case ( 1 )

       call FILE_EXTERNAL_INPUT_get_dims1D( dim1_size, dim1_max, dim1_S, & ! [OUT]
                                            varname, axistype            ) ! [IN]


       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate &
            .or. ( .not. FILE_EXTERNAL_INPUT_item(nid)%file ) ) then
          FILE_EXTERNAL_INPUT_item(nid)%var_size (1) = dim1_size
          FILE_EXTERNAL_INPUT_item(nid)%var_start(1) = dim1_S
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (1) = dim1_max
       else
          if ( dim1_max /= dim_size(1) ) then
             LOG_ERROR("FILE_EXTERNAL_INPUT_regist_var",*) 'data length does not match! ', trim(axistype), ' item:', trim(varname)
             LOG_ERROR_CONT(*) 'dim 1 (data,requested)    : ', dim_size(1), dim1_max
             call PRC_abort
          endif
          FILE_EXTERNAL_INPUT_item(nid)%var_size (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%var_start(1) = 1
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (1) = dim1_max
       end if

       FILE_EXTERNAL_INPUT_item(nid)%transpose    = .false.
       FILE_EXTERNAL_INPUT_item(nid)%dim_size (1) = dim1_size
       FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) = dim1_S
       FILE_EXTERNAL_INPUT_item(nid)%dim_max  (1) = dim1_max

    case ( 2 )

       call FILE_EXTERNAL_INPUT_get_dims2D( dim1_size, dim1_max, dim1_S,             & ! [OUT]
                                            dim2_size, dim2_max, dim2_S,             & ! [OUT]
                                            FILE_EXTERNAL_INPUT_item(nid)%transpose, & ! [OUT]
                                            varname, axistype                        ) ! [IN]

       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate &
            .or. ( .not. FILE_EXTERNAL_INPUT_item(nid)%file ) ) then
          FILE_EXTERNAL_INPUT_item(nid)%var_size (1) = dim1_size
          FILE_EXTERNAL_INPUT_item(nid)%var_start(1) = dim1_S
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%var_size (2) = dim2_size
          FILE_EXTERNAL_INPUT_item(nid)%var_start(2) = dim2_S
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (2) = dim2_max
       else
          if ( dim1_max /= dim_size(1) .OR. dim2_max /= dim_size(2) ) then
             LOG_ERROR("FILE_EXTERNAL_INPUT_regist_var",*) 'data length does not match! ', trim(axistype), ' item:', trim(varname)
             LOG_ERROR_CONT(*) 'dim 1 (data,requested)    : ', dim_size(1), dim1_max
             LOG_ERROR_CONT(*) 'dim 2 (data,requested)    : ', dim_size(2), dim2_max
             call PRC_abort
          endif
          FILE_EXTERNAL_INPUT_item(nid)%var_size (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%var_start(1) = 1
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%var_size (2) = dim2_max
          FILE_EXTERNAL_INPUT_item(nid)%var_start(2) = 1
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (2) = dim2_max
       end if

       if ( FILE_EXTERNAL_INPUT_item(nid)%transpose ) then
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (1) = dim2_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) = dim2_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (1) = dim2_max
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (2) = dim1_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) = dim1_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (2) = dim1_max
       else
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (1) = dim1_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) = dim1_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (2) = dim2_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) = dim2_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (2) = dim2_max
       end if


    case ( 3 )

       call FILE_EXTERNAL_INPUT_get_dims3D( dim1_size, dim1_max, dim1_S,             & ! [OUT]
                                            dim2_size, dim2_max, dim2_S,             & ! [OUT]
                                            dim3_size, dim3_max, dim3_S,             & ! [OUT]
                                            FILE_EXTERNAL_INPUT_item(nid)%transpose, & ! [OUT]
                                            varname, axistype                        ) ! [IN]

       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate &
            .or. ( .not. FILE_EXTERNAL_INPUT_item(nid)%file ) ) then
          FILE_EXTERNAL_INPUT_item(nid)%var_size (1) = dim1_size
          FILE_EXTERNAL_INPUT_item(nid)%var_start(1) = dim1_S
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%var_size (2) = dim2_size
          FILE_EXTERNAL_INPUT_item(nid)%var_start(2) = dim2_S
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (2) = dim2_max
          FILE_EXTERNAL_INPUT_item(nid)%var_size (3) = dim3_size
          FILE_EXTERNAL_INPUT_item(nid)%var_start(3) = dim3_S
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (3) = dim3_max
       else
          if ( dim1_max /= dim_size(1) .OR. dim2_max /= dim_size(2) .OR. dim3_max /= dim_size(3) ) then
             LOG_ERROR("FILE_EXTERNAL_INPUT_regist_var",*) 'data length does not match! ', trim(axistype), ' item:', trim(varname)
             LOG_ERROR_CONT(*) 'dim 1 (data,requested)    : ', dim_size(1), dim1_max
             LOG_ERROR_CONT(*) 'dim 2 (data,requested)    : ', dim_size(2), dim2_max
             LOG_ERROR_CONT(*) 'dim 3 (data,requested)    : ', dim_size(3), dim3_max
             call PRC_abort
          endif
          FILE_EXTERNAL_INPUT_item(nid)%var_size (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%var_start(1) = 1
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%var_size (2) = dim2_max
          FILE_EXTERNAL_INPUT_item(nid)%var_start(2) = 1
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (2) = dim2_max
          FILE_EXTERNAL_INPUT_item(nid)%var_size (3) = dim3_max
          FILE_EXTERNAL_INPUT_item(nid)%var_start(3) = 1
          FILE_EXTERNAL_INPUT_item(nid)%var_max  (3) = dim3_max
       end if

       if ( FILE_EXTERNAL_INPUT_item(nid)%transpose ) then
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (1) = dim3_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) = dim3_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (1) = dim3_max
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (2) = dim1_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) = dim1_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (2) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (3) = dim2_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) = dim2_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (3) = dim2_max
       else
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (1) = dim1_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) = dim1_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (1) = dim1_max
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (2) = dim2_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) = dim2_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (2) = dim2_max
          FILE_EXTERNAL_INPUT_item(nid)%dim_size (3) = dim3_size
          FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) = dim3_S
          FILE_EXTERNAL_INPUT_item(nid)%dim_max  (3) = dim3_max
       end if

    case default
       LOG_ERROR("FILE_EXTERNAL_INPUT_regist_var",*) 'Unexpected dim rank: ', dim_rank
       call PRC_abort
    end select

    do n = dim_rank+1, 3
       FILE_EXTERNAL_INPUT_item(nid)%dim_size (n) = 1
       FILE_EXTERNAL_INPUT_item(nid)%dim_start(n) = 1
       FILE_EXTERNAL_INPUT_item(nid)%dim_max  (n) = 0
       FILE_EXTERNAL_INPUT_item(nid)%var_size (n) = 1
       FILE_EXTERNAL_INPUT_item(nid)%var_start(n) = 1
       FILE_EXTERNAL_INPUT_item(nid)%var_max  (n) = 0
    enddo

    FILE_EXTERNAL_INPUT_item(nid)%ndim        = dim_rank
    FILE_EXTERNAL_INPUT_item(nid)%step_num    = step_num


    allocate( FILE_EXTERNAL_INPUT_item(nid)%value(FILE_EXTERNAL_INPUT_item(nid)%dim_size(1),FILE_EXTERNAL_INPUT_item(nid)%dim_size(2),FILE_EXTERNAL_INPUT_item(nid)%dim_size(3),2) )
    FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,:) = defval

    allocate( FILE_EXTERNAL_INPUT_item(nid)%time(step_num*file_num) )
    do n = 1, FILE_EXTERNAL_INPUT_item(nid)%step_num*file_num
       FILE_EXTERNAL_INPUT_item(nid)%time(n) = time_step * ( n - 1 ) + time_now
    end do

    FILE_EXTERNAL_INPUT_item(nid)%data_step_prev = 1
    FILE_EXTERNAL_INPUT_item(nid)%data_step_next = 2

    return
  end subroutine FILE_EXTERNAL_INPUT_regist_var

  !-----------------------------------------------------------------------------
  !> Regist external data
  subroutine FILE_EXTERNAL_INPUT_regist_external_1D( &
       varname,      &
       var,          &
       axistype,     &
       step_nmax,    &
       time_now,     &
       time_step,    &
       aggregate,    &
       allow_missing )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    character(len=*), intent(in) :: varname
    real(RP),         intent(in) :: var(:)
    character(len=*), intent(in) :: axistype
    integer,          intent(in) :: step_nmax
    real(DP),         intent(in) :: time_now
    real(DP),         intent(in) :: time_step
    logical,          intent(in),  optional :: aggregate
    logical,          intent(in),  optional :: allow_missing

    integer :: dim_size(1)
    logical :: error
    !---------------------------------------------------------------------------

    call FILE_EXTERNAL_INPUT_init_var( &
         varname,      &
         axistype,     &
         aggregate,    &
         allow_missing )

    dim_size(:) = shape( var )

    call FILE_EXTERNAL_INPUT_regist_var( &
         varname,     &
         axistype,    &
         1,           &
         dim_size(:), &
         step_nmax,   &
         1,           &
         time_now,    &
         time_step,   &
         UNDEF        )

    call FILE_EXTERNAL_INPUT_put_ref_1D( &
         varname, &
         var(:),  &
         error    )

    return
  end subroutine FILE_EXTERNAL_INPUT_regist_external_1D

  !-----------------------------------------------------------------------------
  !> Regist external data
  subroutine FILE_EXTERNAL_INPUT_regist_external_2D( &
       varname,      &
       var,     &
       axistype,     &
       step_nmax,    &
       time_now,     &
       time_step,    &
       aggregate,    &
       allow_missing )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    character(len=*), intent(in) :: varname
    real(RP),         intent(in) :: var(:,:)
    character(len=*), intent(in) :: axistype
    integer,          intent(in) :: step_nmax
    real(DP),         intent(in) :: time_now
    real(DP),         intent(in) :: time_step
    logical,          intent(in),  optional :: aggregate
    logical,          intent(in),  optional :: allow_missing

    integer :: dim_size(2)
    logical :: error
    !---------------------------------------------------------------------------

    call FILE_EXTERNAL_INPUT_init_var( &
         varname,      &
         axistype,     &
         aggregate,    &
         allow_missing )

    dim_size(:) = shape( var )

    call FILE_EXTERNAL_INPUT_regist_var( &
         varname,     &
         axistype,    &
         2,           &
         dim_size(:), &
         step_nmax,   &
         1,           &
         time_now,    &
         time_step,   &
         UNDEF        )

    call FILE_EXTERNAL_INPUT_put_ref_2D( &
         varname,  &
         var(:,:), &
         error     )

    return
  end subroutine FILE_EXTERNAL_INPUT_regist_external_2D

  !-----------------------------------------------------------------------------
  !> Regist external data
  subroutine FILE_EXTERNAL_INPUT_regist_external_3D( &
       varname,      &
       var,          &
       axistype,     &
       step_nmax,    &
       time_now,     &
       time_step,    &
       aggregate,    &
       allow_missing )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    character(len=*), intent(in) :: varname
    real(RP),         intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: axistype
    integer,          intent(in) :: step_nmax
    real(DP),         intent(in) :: time_now
    real(DP),         intent(in) :: time_step
    logical,          intent(in),  optional :: aggregate
    logical,          intent(in),  optional :: allow_missing

    integer :: dim_size(3)
    logical :: error
    !---------------------------------------------------------------------------

    call FILE_EXTERNAL_INPUT_init_var( &
         varname,      &
         axistype,     &
         aggregate,    &
         allow_missing )

    dim_size(:) = shape( var )

    call FILE_EXTERNAL_INPUT_regist_var( &
         varname,     &
         axistype,    &
         3,           &
         dim_size(:), &
         step_nmax,   &
         1,           &
         time_now,    &
         time_step,   &
         UNDEF        )

    call FILE_EXTERNAL_INPUT_put_ref_3D( &
         varname,    &
         var(:,:,:), &
         error       )

    return
  end subroutine FILE_EXTERNAL_INPUT_regist_external_3D

  !-----------------------------------------------------------------------------
  !> Regist data
  subroutine FILE_EXTERNAL_INPUT_regist_file( &
       basename,              &
       basename_add_num,      &
       number_of_files,       &
       varname,               &
       axistype,              &
       enable_periodic_year,  &
       enable_periodic_month, &
       enable_periodic_day,   &
       step_fixed,            &
       defval,                &
       check_coordinates,     &
       aggregate,             &
       allow_missing,         &
       step_limit,            &
       update_dt,             &
       exist                  )
    use scale_file_h, only: &
       FILE_FREAD
    use scale_file, only: &
       FILE_AGGREGATE,        &
       FILE_Open,             &
       FILE_Get_All_DataInfo, &
       FILE_Read
    use scale_prc, only: &
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

    character(len=*), intent(in)  :: basename
    logical,          intent(in)  :: basename_add_num
    integer,          intent(in)  :: number_of_files
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype
    integer,          intent(in)  :: step_fixed            ! fixed step position to read
    logical,          intent(in)  :: enable_periodic_year  ! treat as yearly               periodic data?
    logical,          intent(in)  :: enable_periodic_month ! treat as yearly,monthly       periodic data?
    logical,          intent(in)  :: enable_periodic_day   ! treat as yearly,monthly,daily periodic data?
    real(RP),         intent(in)  :: defval

    logical,          intent(in),  optional :: check_coordinates
    logical,          intent(in),  optional :: aggregate
    logical,          intent(in),  optional :: allow_missing
    integer,          intent(in),  optional :: step_limit            ! limit number for reading data
    real(DP),         intent(out), optional :: update_dt
    logical,          intent(out), optional :: exist

    integer                :: step_nmax
    character(len=H_MID)   :: description
    character(len=H_SHORT) :: unit
    character(len=H_MID)   :: standard_name
    integer                :: datatype
    integer                :: dim_rank
    character(len=H_SHORT) :: dim_name  (FILE_EXTERNAL_INPUT_dim_limit)
    integer                :: dim_size  (FILE_EXTERNAL_INPUT_dim_limit)
    integer                :: var_size  (FILE_EXTERNAL_INPUT_dim_limit)
    integer                :: natts
    character(len=H_SHORT) :: att_name  (FILE_EXTERNAL_INPUT_att_limit)
    integer                :: att_type  (FILE_EXTERNAL_INPUT_att_limit)
    integer                :: att_len   (FILE_EXTERNAL_INPUT_att_limit)
    real(DP)               :: time_start(FILE_EXTERNAL_INPUT_step_limit)
    real(DP)               :: time_end  (FILE_EXTERNAL_INPUT_step_limit)
    character(len=H_MID)   :: time_units
    character(len=H_SHORT) :: calendar

    integer  :: datadate(6)   !< date
    real(DP) :: datasubsec    !< subsecond
    integer  :: dataday       !< absolute day
    real(DP) :: datasec       !< absolute second
    integer  :: offset_year   !< offset year

    character(len=H_LONG) :: filename

    logical :: aggregate_
    integer :: step_limit_

    real(RP), allocatable :: buf(:,:,:)

    logical :: error

    integer :: fid
    integer :: nid, n
    !---------------------------------------------------------------------------

    if ( present(aggregate) ) then
       aggregate_ = aggregate
    else
       aggregate_ = FILE_AGGREGATE
    end if

    if ( present(step_limit) ) then
       if ( step_limit > 0 ) then
          step_limit_ = step_limit
       else
          step_limit_ = FILE_EXTERNAL_INPUT_step_limit
       endif
    else
       step_limit_ = FILE_EXTERNAL_INPUT_step_limit
    endif

    call FILE_EXTERNAL_INPUT_init_var( &
         varname,      &
         axistype,     &
         aggregate_,   &
         allow_missing )


    nid = FILE_EXTERNAL_INPUT_item_count

    FILE_EXTERNAL_INPUT_item(nid)%file             = .true.
    FILE_EXTERNAL_INPUT_item(nid)%nfile            = number_of_files
    FILE_EXTERNAL_INPUT_item(nid)%file_current     = 1
    FILE_EXTERNAL_INPUT_item(nid)%data_step_offset = 0
    FILE_EXTERNAL_INPUT_item(nid)%step_limit       = step_limit_

    allocate( FILE_EXTERNAL_INPUT_item(nid)%basename(number_of_files) )
    if ( number_of_files > 1 .or. basename_add_num ) then
       do n = 1, number_of_files
          write(filename,'(A,A,I5.5)') trim(basename), '_', n - 1
          FILE_EXTERNAL_INPUT_item(nid)%basename(n) = filename
       enddo
    else
       FILE_EXTERNAL_INPUT_item(nid)%basename(1) = basename
    end if


    filename   = FILE_EXTERNAL_INPUT_item(nid)%basename(1)
    call FILE_Open( filename,             & ! [IN]
                    fid,                  & ! [OUT]
                    aggregate=aggregate_, & ! [IN]
                    rankid=PRC_myrank     ) ! [IN]

    ! read from file
    call FILE_Get_All_Datainfo( fid, varname,                                       & ! [IN]
                                step_nmax,                                          & ! [OUT]
                                description, unit,  standard_name,                  & ! [OUT]
                                datatype,                                           & ! [OUT]
                                dim_rank, dim_name(:), dim_size(:),                 & ! [OUT]
                                natts, att_name(:), att_type(:), att_len(:),        & ! [OUT]
                                time_start(1:step_limit_), time_end(1:step_limit_), & ! [OUT]
                                time_units, calendar                                ) ! [OUT]

    if ( step_nmax > 0 ) then
       if ( present(exist) ) then
          exist = .true.
       endif
    else
       if ( present(exist) ) then
          exist = .false.
          return
       else
          LOG_ERROR("FILE_EXTERNAL_INPUT_regist",*) 'Data not found! filename,varname = ', trim(filename), ', ', trim(varname)
          call PRC_abort
       endif
    endif

    do n = dim_rank+1, 3
       dim_size(n) = 1
    enddo

    FILE_EXTERNAL_INPUT_item(nid)%fid           = fid

    if ( enable_periodic_day ) then
       FILE_EXTERNAL_INPUT_item(nid)%flag_periodic = I_periodic_day
    elseif( enable_periodic_month ) then
       FILE_EXTERNAL_INPUT_item(nid)%flag_periodic = I_periodic_month
    elseif( enable_periodic_year ) then
       FILE_EXTERNAL_INPUT_item(nid)%flag_periodic = I_periodic_year
    endif

    call FILE_EXTERNAL_INPUT_regist_var( &
         varname,         &
         axistype,        &
         dim_rank,        &
         dim_size(:),     &
         step_nmax,       &
         number_of_files, &
         0.0_DP, 0.0_DP,  &
         defval           )

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

             LOG_INFO("FILE_EXTERNAL_INPUT_regist",*) 'data time is updated.'
          endif

       else ! normal mode

          if (      FILE_EXTERNAL_INPUT_item(nid)%data_step_next == 1                          &
               .OR. FILE_EXTERNAL_INPUT_item(nid)%data_step_next == FILE_EXTERNAL_INPUT_item(nid)%step_num+1 ) then
             LOG_ERROR("FILE_EXTERNAL_INPUT_regist",*) 'Current time is out of period of external data! ', trim(varname)
             call PRC_abort
          endif

       endif

    endif

    !--- read first data
    LOG_INFO("FILE_EXTERNAL_INPUT_regist",'(1x,A,A15)') 'Initial read of external data : ', trim(varname)

    allocate( buf(FILE_EXTERNAL_INPUT_item(nid)%var_size(1),FILE_EXTERNAL_INPUT_item(nid)%var_size(2),FILE_EXTERNAL_INPUT_item(nid)%var_size(3)) )

    select case ( dim_rank )
    case ( 1 )

       ! read prev
       LOG_INFO("FILE_EXTERNAL_INPUT_regist",'(1x,A,A,A,I4,A)') &
                  'Read 1D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_prev, ')'

       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
          call FILE_EXTERNAL_INPUT_read_1d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                            buf(:,1,1),                                       & ! [OUT]
                                            step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev ) ! [IN]
       else
          call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                          FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                          buf(:,1,1),                                       & ! [OUT]
                          step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev ) ! [IN]
       end if

       call FILE_EXTERNAL_INPUT_put_ref_1D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                            buf(:,1,1),                            &
                                            error                                  )

       ! read next
       LOG_INFO("FILE_EXTERNAL_INPUT_regist",'(1x,A,A,A,I4,A)') &
                'Read 1D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ')'

       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
          call FILE_EXTERNAL_INPUT_read_1d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                            buf(:,1,1),                                       & ! [OUT]
                                            step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
       else
          call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                          FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                          buf(:,1,1),                                       & ! [OUT]
                          step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
       end if

       call FILE_EXTERNAL_INPUT_put_ref_1D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                            buf(:,1,1),                            &
                                            error                                  )

    case ( 2 )

       ! read prev
       LOG_INFO("FILE_EXTERNAL_INPUT_regist",'(1x,A,A,A,I4,A)') &
                  'Read 2D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_prev, ')'

       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
          call FILE_EXTERNAL_INPUT_read_2d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                            buf(:,:,1),                                       & ! [OUT]
                                            step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev ) ! [IN]
       else
          call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                          FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                          buf(:,:,1),                                       & ! [OUT]
                          step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev ) ! [IN]
       end if

       call FILE_EXTERNAL_INPUT_put_ref_2D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                            buf(:,:,1),                            &
                                            error                                  )

       ! read next
       LOG_INFO("FILE_EXTERNAL_INPUT_regist",'(1x,A,A,A,I4,A)') &
                  'Read 2D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ')'

       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
          call FILE_EXTERNAL_INPUT_read_2d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                            buf(:,:,1),                                       & ! [OUT]
                                            step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
       else
          call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                          FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                          buf(:,:,1),                                       & ! [OUT]
                          step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
       end if

       call FILE_EXTERNAL_INPUT_put_ref_2D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                            buf(:,:,1),                            &
                                            error                                  )

    case ( 3 )

       ! read prev
       LOG_INFO("FILE_EXTERNAL_INPUT_regist",'(1x,A,A,A,I4,A)') &
                  'Read 3D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_prev, ')'

       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
          call FILE_EXTERNAL_INPUT_read_3d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                            buf(:,:,:),                                       & ! [OUT]
                                            step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev ) ! [IN]
       else
          call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                          FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                          buf(:,:,:),                                       & ! [OUT]
                          step=FILE_EXTERNAL_INPUT_item(nid)%data_step_prev ) ! [IN]
       end if

       call FILE_EXTERNAL_INPUT_put_ref_3D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                            buf(:,:,:),                            &
                                            error                                  )

       ! read next
       LOG_INFO("FILE_EXTERNAL_INPUT_regist",'(1x,A,A,A,I4,A)') &
                  'Read 3D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
                  ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ')'

       if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
          call FILE_EXTERNAL_INPUT_read_3d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                            FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                            buf(:,:,:),                                       & ! [OUT]
                                            step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
       else
          call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                          FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                          buf(:,:,:),                                       & ! [OUT]
                          step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
       end if

       call FILE_EXTERNAL_INPUT_put_ref_3D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                            buf(:,:,:),                            &
                                            error                                  )

    case default
       LOG_ERROR("FILE_EXTERNAL_INPUT_regist",*) 'Unexpected dim rank: ', dim_rank
       call PRC_abort
    end select

    deallocate( buf )

    if ( present(check_coordinates) ) then
       if ( check_coordinates ) then
          call FILE_CARTESC_check_coordinates( fid,                                                &
                                               atmos     = FILE_EXTERNAL_INPUT_item(nid)%ndim==3,  &
                                               transpose = FILE_EXTERNAL_INPUT_item(nid)%transpose )
       endif
    endif


    if ( present(update_dt) ) then
       update_dt = FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_next) &
                 - FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_prev)
    end if

    return
  end subroutine FILE_EXTERNAL_INPUT_regist_file

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine FILE_EXTERNAL_INPUT_update_1D( &
       varname,      &
       time_current, &
       var,          &
       error         )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_file, only: &
       FILE_Read
    implicit none
    character(len=*), intent(in)  :: varname      ! item name
    real(DP),         intent(in)  :: time_current ! current time
    real(RP),         intent(out) :: var(:)       ! variable
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile
    integer  :: step_next

    real(RP), allocatable :: buf(:)

    integer  :: n
    integer  :: n1, nn1
    !---------------------------------------------------------------------------

    nid = FILE_EXTERNAL_INPUT_getId(varname)

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_update_1D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    if ( FILE_EXTERNAL_INPUT_item(nid)%ndim /= 1 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_update_1D",*) 'Data is not 1D var: ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
       error = .true.
       return
    endif

    call FILE_EXTERNAL_INPUT_time_advance( nid,          & ! [IN]
                                           time_current, & ! [IN]
                                           weight,       & ! [OUT]
                                           do_readfile   ) ! [OUT]

    if ( do_readfile ) then

       if ( FILE_EXTERNAL_INPUT_item(nid)%file ) then

          ! next -> prev
          FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

          step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next - FILE_EXTERNAL_INPUT_item(nid)%data_step_offset

          LOG_INFO("FILE_EXTERNAL_INPUT_update_1D",'(1x,A,A,A,I4,A,I4,A)') &
               'Read 1D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
               ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ', file step=', step_next, ')'

          allocate( buf(FILE_EXTERNAL_INPUT_item(nid)%var_size(1)) )

          ! read next
          if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
             call FILE_EXTERNAL_INPUT_read_1d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                               FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                               FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                               buf(:),                                           & ! [OUT]
                                               step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
          else
             call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,     & ! [IN]
                             FILE_EXTERNAL_INPUT_item(nid)%varname, & ! [IN]
                             buf(:),                                & ! [OUT]
                             step=step_next                         ) ! [IN]
          end if

          call FILE_EXTERNAL_INPUT_put_ref_1D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                               buf(:),                                &
                                               error                                  )

          deallocate( buf )

       end if

    endif


    error = .false.

    ! store data with weight
    do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(1)
       nn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1

       if (       abs( FILE_EXTERNAL_INPUT_item(nid)%value(nn1,1,1,I_prev) - UNDEF ) > abs( UNDEF * 0.1_RP ) &
            .and. abs( FILE_EXTERNAL_INPUT_item(nid)%value(nn1,1,1,I_next) - UNDEF ) > abs( UNDEF * 0.1_RP ) ) then
          var(nn1) = ( 1.0_RP-weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(nn1,1,1,I_prev) &
                   + (        weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(nn1,1,1,I_next)
       else
          if ( FILE_EXTERNAL_INPUT_item(nid)%allow_missing ) then
             var(nn1) = UNDEF
          else
             LOG_INFO("FILE_EXTERNAL_INPUT_update_1D",*) 'missing value is found in ', &
                  trim(FILE_EXTERNAL_INPUT_item(nid)%varname), ' at (',nn1,')'
             error = .true.
             exit
          end if
       end if
    enddo

    return
  end subroutine FILE_EXTERNAL_INPUT_update_1D

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine FILE_EXTERNAL_INPUT_update_2D( &
       varname,      &
       time_current, &
       var,          &
       error         )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
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

    real(RP), allocatable :: buf(:,:)

    integer :: n
    integer :: n1, n2
    integer :: nn1, nn2
    !---------------------------------------------------------------------------

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_update_2D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    if ( FILE_EXTERNAL_INPUT_item(nid)%ndim /= 2 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_update_2D",*) 'Data is not 2D var: ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
       error = .true.
       return
    endif

    call FILE_EXTERNAL_INPUT_time_advance( nid,          & ! [IN]
                                           time_current, & ! [IN]
                                           weight,       & ! [OUT]
                                           do_readfile   ) ! [OUT]

    if ( do_readfile ) then

       if ( FILE_EXTERNAL_INPUT_item(nid)%file ) then

          ! next -> prev
          FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

          step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next - FILE_EXTERNAL_INPUT_item(nid)%data_step_offset

          LOG_INFO("FILE_EXTERNAL_INPUT_update_2D",'(1x,A,A,A,I4,A,I4,A)') &
               'Read 2D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
               ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ', file step=', step_next, ')'

          allocate( buf(FILE_EXTERNAL_INPUT_item(nid)%var_size(1),FILE_EXTERNAL_INPUT_item(nid)%var_size(2)) )

          ! read next
          if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
             call FILE_EXTERNAL_INPUT_read_2d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                               FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                               FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                               buf(:,:),                                         & ! [OUT]
                                               step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
          else
             call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,     & ! [IN]
                             FILE_EXTERNAL_INPUT_item(nid)%varname, & ! [IN]
                             buf(:,:),                              & ! [OUT]
                             step=step_next                         ) ! [IN]
          end if

          call FILE_EXTERNAL_INPUT_put_ref_2D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                               buf(:,:),                              &
                                               error                                  )

          deallocate( buf )

       end if

    endif

    error = .false.

    ! store data with weight
    do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(2)
       nn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1

       do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(1)
          nn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1

          if (       abs( FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,1,I_prev) - UNDEF ) > abs( UNDEF * 0.1_RP ) &
               .and. abs( FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,1,I_next) - UNDEF ) > abs( UNDEF * 0.1_RP ) ) then
             var(nn1,nn2) = ( 1.0_RP-weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,1,I_prev) &
                          + (        weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,1,I_next)
          else
             if ( FILE_EXTERNAL_INPUT_item(nid)%allow_missing ) then
                var(nn1,nn2) = UNDEF
             else
                LOG_INFO("FILE_EXTERNAL_INPUT_update_2D",*) 'missing value is found in ', &
                     trim(FILE_EXTERNAL_INPUT_item(nid)%varname), ' at (',nn1,',',nn2,')'
                error = .true.
                exit
             end if
          end if
       enddo
    enddo

    return
  end subroutine FILE_EXTERNAL_INPUT_update_2D

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine FILE_EXTERNAL_INPUT_update_3D( &
       varname,      &
       time_current, &
       var,          &
       error         )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
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

    real(RP), allocatable :: buf(:,:,:)

    integer :: n
    integer :: n1, n2, n3
    integer :: nn1, nn2, nn3
    !---------------------------------------------------------------------------

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_update_3D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    if ( FILE_EXTERNAL_INPUT_item(nid)%ndim /= 3 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_update_3D",*) 'Data is not 3D var: ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
       error = .true.
       return
    endif

    call FILE_EXTERNAL_INPUT_time_advance( nid,          & ! [IN]
                                           time_current, & ! [IN]
                                           weight,       & ! [OUT]
                                           do_readfile   ) ! [OUT]

    if ( do_readfile ) then

       if ( FILE_EXTERNAL_INPUT_item(nid)%file ) then

          ! next -> prev
          FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

          step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next - FILE_EXTERNAL_INPUT_item(nid)%data_step_offset

          LOG_INFO("FILE_EXTERNAL_INPUT_update_3D",'(1x,A,A,A,I4,A,I4,A)') &
               'Read 3D var           : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname), &
               ' (step= ', FILE_EXTERNAL_INPUT_item(nid)%data_step_next, ', file step=', step_next, ')'

          allocate( buf(FILE_EXTERNAL_INPUT_item(nid)%var_size(1),FILE_EXTERNAL_INPUT_item(nid)%var_size(2),FILE_EXTERNAL_INPUT_item(nid)%var_size(3)) )

          ! read next
          if ( FILE_EXTERNAL_INPUT_item(nid)%aggregate ) then
             call FILE_EXTERNAL_INPUT_read_3d( FILE_EXTERNAL_INPUT_item(nid)%fid,                & ! [IN]
                                               FILE_EXTERNAL_INPUT_item(nid)%varname,            & ! [IN]
                                               FILE_EXTERNAL_INPUT_item(nid)%axistype,           & ! [IN]
                                               buf(:,:,:),                                       & ! [OUT]
                                               step=FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ! [IN]
          else
             call FILE_Read( FILE_EXTERNAL_INPUT_item(nid)%fid,     & ! [IN]
                             FILE_EXTERNAL_INPUT_item(nid)%varname, & ! [IN]
                             buf(:,:,:),                            & ! [OUT]
                             step=step_next                         ) ! [IN]
          end if

          call FILE_EXTERNAL_INPUT_put_ref_3D( FILE_EXTERNAL_INPUT_item(nid)%varname, &
                                               buf(:,:,:),                            &
                                               error                                  )

       end if

    endif

    error = .false.

    ! store data with weight
    do n3 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(3)
       nn3 = n3 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) - 1

       do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(2)
          nn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1

          do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(1)
             nn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1

             if (       abs( FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,nn3,I_prev) - UNDEF ) > abs( UNDEF * 0.1_RP ) &
                  .and. abs( FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,nn3,I_next) - UNDEF ) > abs( UNDEF * 0.1_RP ) ) then
                var(nn1,nn2,nn3) = ( 1.0_RP-weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,nn3,I_prev) &
                                 + (        weight ) * FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,nn3,I_next)
             else
                if ( FILE_EXTERNAL_INPUT_item(nid)%allow_missing ) then
                   var(nn1,nn2,nn3) = UNDEF
                else
                   LOG_INFO("FILE_EXTERNAL_INPUT_update_3D",*) 'missing value is found in ', &
                        trim(FILE_EXTERNAL_INPUT_item(nid)%varname), ' at (',nn1,',',nn2,',',nn3,')'
                   error = .true.
                   exit
                end if
             end if
          enddo
       enddo
    enddo

    return
  end subroutine FILE_EXTERNAL_INPUT_update_3D

  !-----------------------------------------------------------------------------
  !> Put reference data
  subroutine FILE_EXTERNAL_INPUT_put_ref_1D( &
       varname, &
       var,     &
       error    )
    implicit none
    character(len=*), intent(in)  :: varname ! item name
    real(RP),         intent(in)  :: var(:)  ! variable
    logical,          intent(out) :: error   ! error code

    integer :: nid
    integer :: n1, nn1, nnn1

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_put_ref_1D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    error = .false.

    FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

    do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(1)
       nn1  = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1
       nnn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%var_start(1) - 1
       FILE_EXTERNAL_INPUT_item(nid)%value(nn1,1,1,I_next) = var(nnn1)
    enddo

    return
  end subroutine FILE_EXTERNAL_INPUT_put_ref_1D

  subroutine FILE_EXTERNAL_INPUT_put_ref_2D( &
       varname, &
       var,     &
       error    )
    implicit none
    character(len=*), intent(in)  :: varname  ! item name
    real(RP),         intent(in)  :: var(:,:) ! variable
    logical,          intent(out) :: error    ! error code

    integer :: nid
    integer :: n1, n2, nn1, nn2, nnn1, nnn2

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_put_var_2D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    error = .false.

    FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

    if ( FILE_EXTERNAL_INPUT_item(nid)%transpose ) then
       ! (x,z)->(z,x)
       do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(2)
          nn2  = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1
          nnn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%var_start(1) - 1

          do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(1)
             nn1  = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1
             nnn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%var_start(2) - 1
             FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,1,I_next) = var(nnn2,nnn1)
          enddo
       enddo
    else
       ! (z,x)->(z,x)
       do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(2)
          nn2  = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1
          nnn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%var_start(2) - 1

          do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(1)
             nn1  = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1
             nnn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%var_start(1) - 1
             FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,1,I_next) = var(nnn1,nnn2)
          enddo
       enddo
    endif

    return
  end subroutine FILE_EXTERNAL_INPUT_put_ref_2D

  subroutine FILE_EXTERNAL_INPUT_put_ref_3D( &
       varname, &
       var,     &
       error    )
    implicit none
    character(len=*), intent(in)  :: varname    ! item name
    real(RP),         intent(in)  :: var(:,:,:) ! variable
    logical,          intent(out) :: error      ! error code

    integer :: nid
    integer :: n1, n2, n3, nn1, nn2, nn3, nnn1, nnn2, nnn3

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_put_ref_3D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    error = .false.

    FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_prev) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,I_next)

    if ( FILE_EXTERNAL_INPUT_item(nid)%transpose ) then
       ! (x,y,z)->(z,x,y)
       do n3 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(3)
          nn3  = n3 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) - 1
          nnn3 = n3 + FILE_EXTERNAL_INPUT_item(nid)%var_start(2) - 1

          do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(2)
             nn2  = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1
             nnn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%var_start(1) - 1

             do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(1)
                nn1  = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1
                nnn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%var_start(3) - 1
                FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,nn3,I_next) = var(nnn2,nnn3,nnn1)
             enddo
          enddo
       enddo
    else
       ! (z,x,y)->(z,x,y)
       do n3 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(3)
          nn3  = n3 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(3) - 1
          nnn3 = n3 + FILE_EXTERNAL_INPUT_item(nid)%var_start(3) - 1

          do n2 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(2)
             nn2  = n2 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(2) - 1
             nnn2 = n2 + FILE_EXTERNAL_INPUT_item(nid)%var_start(2) - 1

             do n1 = 1, FILE_EXTERNAL_INPUT_item(nid)%dim_max(1)
                nn1  = n1 + FILE_EXTERNAL_INPUT_item(nid)%dim_start(1) - 1
                nnn1 = n1 + FILE_EXTERNAL_INPUT_item(nid)%var_start(1) - 1
                FILE_EXTERNAL_INPUT_item(nid)%value(nn1,nn2,nn3,I_next) = var(nnn1,nnn2,nnn3)
             enddo
          enddo
       enddo
    endif

    return
  end subroutine FILE_EXTERNAL_INPUT_put_ref_3D

  !-----------------------------------------------------------------------------
  !> Get reference data
  subroutine FILE_EXTERNAL_INPUT_get_ref_1D( &
       varname, &
       var,     &
       error,   &
       i_step   )
    implicit none
    character(len=*), intent(in)  :: varname ! item name
    real(RP),         intent(out) :: var(:)  ! variable
    logical,          intent(out) :: error   ! error code

    integer, optional, intent(in)  :: i_step

    integer :: i_step_
    integer :: nid

    if ( present(i_step) ) then
       i_step_ = i_step
    else
       i_step_ = I_next
    end if

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_get_ref_3D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    error = .false.

    var(:) = FILE_EXTERNAL_INPUT_item(nid)%value(:,1,1,i_step_)

    return
  end subroutine FILE_EXTERNAL_INPUT_get_ref_1D

  subroutine FILE_EXTERNAL_INPUT_get_ref_2D( &
       varname, &
       var,     &
       error,   &
       i_step   )
    implicit none
    character(len=*), intent(in)  :: varname  ! item name
    real(RP),         intent(out) :: var(:,:) ! variable
    logical,          intent(out) :: error    ! error code

    integer, optional, intent(in)  :: i_step

    integer :: i_step_
    integer :: nid

    if ( present(i_step) ) then
       i_step_ = i_step
    else
       i_step_ = I_next
    end if

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_update_3D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    error = .false.

    var(:,:) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,1,i_step_)

    return
  end subroutine FILE_EXTERNAL_INPUT_get_ref_2D

  !-----------------------------------------------------------------------------
  !> Get reference data
  subroutine FILE_EXTERNAL_INPUT_get_ref_3D( &
       varname, &
       var,     &
       error,   &
       i_step   )
    implicit none
    character(len=*), intent(in)  :: varname    ! item name
    real(RP),         intent(out) :: var(:,:,:) ! variable
    logical,          intent(out) :: error      ! error code

    integer, optional, intent(in)  :: i_step

    integer :: i_step_
    integer :: nid

    if ( present(i_step) ) then
       i_step_ = i_step
    else
       i_step_ = I_next
    end if

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    if ( nid == 0 ) then
       LOG_INFO("FILE_EXTERNAL_INPUT_update_3D",*) 'Variable was not registered: ', trim(varname)
       error = .true.
       return
    endif

    error = .false.

    var(:,:,:) = FILE_EXTERNAL_INPUT_item(nid)%value(:,:,:,i_step_)

    return
  end subroutine FILE_EXTERNAL_INPUT_get_ref_3D

  !-----------------------------------------------------------------------------
  !> Check time to read
  subroutine FILE_EXTERNAL_INPUT_query( &
       varname,          &
       time_current, &
       do_readdata   )
    implicit none
    character(len=*), intent(in)  :: varname      ! variable name
    real(DP),         intent(in)  :: time_current ! current time
    logical,          intent(out) :: do_readdata  ! read new data at this time?

    integer :: nid

    nid = FILE_EXTERNAL_INPUT_getId( varname )

    call FILE_EXTERNAL_INPUT_query_Id( nid,          &
                                       time_current, &
                                       do_readdata   )

    return
  end subroutine FILE_EXTERNAL_INPUT_query

  subroutine FILE_EXTERNAL_INPUT_query_Id( &
       nid,          &
       time_current, &
       do_readdata   )
    implicit none
    integer,  intent(in)  :: nid ! variable id
    real(DP), intent(in)  :: time_current ! current time
    logical,  intent(out) :: do_readdata  ! read new data at this time?

    if (       ( .not. FILE_EXTERNAL_INPUT_item(nid)%fixed_step ) &
         .and. ( time_current > FILE_EXTERNAL_INPUT_item(nid)%time( FILE_EXTERNAL_INPUT_item(nid)%data_step_next ) ) &
         ) then
       do_readdata = .true.
    else
       do_readdata = .false.
    end if

    return
  end subroutine FILE_EXTERNAL_INPUT_query_Id

  !-----------------------------------------------------------------------------
  !> Check time to read and calc weight
  subroutine FILE_EXTERNAL_INPUT_time_advance( &
       nid,          &
       time_current, &
       weight,       &
       do_readdata   )
    use scale_file_h, only: &
       FILE_FREAD
    use scale_file, only: &
       FILE_Open,           &
       FILE_Get_All_DataInfo
    use scale_prc, only: &
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
    logical,  intent(out) :: do_readdata  ! read new data at this time?

    integer                :: step_nmax
    character(len=H_MID)   :: description
    character(len=H_SHORT) :: unit
    character(len=H_MID)   :: standard_name
    integer                :: datatype
    integer                :: dim_rank
    character(len=H_SHORT) :: dim_name  (FILE_EXTERNAL_INPUT_dim_limit)
    integer                :: dim_size  (FILE_EXTERNAL_INPUT_dim_limit)
    integer                :: natts
    character(len=H_SHORT) :: att_name  (FILE_EXTERNAL_INPUT_att_limit)
    integer                :: att_type  (FILE_EXTERNAL_INPUT_att_limit)
    integer                :: att_len   (FILE_EXTERNAL_INPUT_att_limit)
    real(DP)               :: time_start(FILE_EXTERNAL_INPUT_step_limit)
    real(DP)               :: time_end  (FILE_EXTERNAL_INPUT_step_limit)
    character(len=H_MID)   :: time_units
    character(len=H_SHORT) :: calendar

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

    call FILE_EXTERNAL_INPUT_query_Id( nid,          & ! [IN]
                                       time_current, & ! [IN]
                                       do_readdata   ) ! [OUT]

    if ( do_readdata ) then

       LOG_INFO("FILE_EXTERNAL_INPUT_time_advance",'(1x,A,A15)') 'Update external input : ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)

       ! update step position
       FILE_EXTERNAL_INPUT_item(nid)%data_step_prev = FILE_EXTERNAL_INPUT_item(nid)%data_step_next
       FILE_EXTERNAL_INPUT_item(nid)%data_step_next = FILE_EXTERNAL_INPUT_item(nid)%data_step_next + 1

       if ( FILE_EXTERNAL_INPUT_item(nid)%file ) then

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
                   call FILE_Get_All_Datainfo( fid, FILE_EXTERNAL_INPUT_item(nid)%varname,             & ! [IN]
                                               step_nmax,                                              & ! [OUT]
                                               description, unit, standard_name,                       & ! [OUT]
                                               datatype,                                               & ! [OUT]
                                               dim_rank, dim_name(:), dim_size(:),                     & ! [OUT]
                                               natts, att_name(:), att_type(:), att_len(:),            & ! [OUT]
                                               time_start(1:FILE_EXTERNAL_INPUT_item(nid)%step_limit), & ! [OUT]
                                               time_end  (1:FILE_EXTERNAL_INPUT_item(nid)%step_limit), & ! [OUT]
                                               time_units, calendar                                    ) ! [OUT]

                   if ( step_nmax == 0 ) then
                      LOG_ERROR("FILE_EXTERNAL_INPUT_time_advance",*) 'Data not found! basename = ', trim(FILE_EXTERNAL_INPUT_item(nid)%basename(FILE_EXTERNAL_INPUT_item(nid)%file_current)), &
                                                    ', varname = ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
                      call PRC_abort
                   endif

                   do n = 1, dim_rank
                      if ( FILE_EXTERNAL_INPUT_item(nid)%var_size(n) /= dim_size(n) ) then
                         LOG_ERROR("FILE_EXTERNAL_INPUT_time_advance",*) 'The size of dimension', n, ' is inconsistent! '
                         LOG_ERROR_CONT(*) 'size (previous,current) = ', FILE_EXTERNAL_INPUT_item(nid)%var_size(n), dim_size(n)
                         LOG_ERROR_CONT(*) 'basename = ', trim(FILE_EXTERNAL_INPUT_item(nid)%basename(FILE_EXTERNAL_INPUT_item(nid)%file_current)), &
                                       ', varname = ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
                         call PRC_abort
                      endif
                   enddo

                   do n = 1, step_nmax
                      nn = FILE_EXTERNAL_INPUT_item(nid)%step_num + n
                      FILE_EXTERNAL_INPUT_item(nid)%time(nn) = CALENDAR_CFunits2sec( time_end(n), time_units, TIME_OFFSET_year, TIME_STARTDAYSEC )
                   enddo

                   if ( FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_prev) > FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_next) ) then
                      LOG_ERROR("FILE_EXTERNAL_INPUT_time_advance",*) 'Time in new file is earlier than last time of previous file! stop'
                      LOG_ERROR_CONT(*) 'Time (previous,current)  = ', FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_prev), &
                                                                    FILE_EXTERNAL_INPUT_item(nid)%time(FILE_EXTERNAL_INPUT_item(nid)%data_step_next)
                      LOG_ERROR_CONT(*) 'Data not found! basename = ', trim(FILE_EXTERNAL_INPUT_item(nid)%basename(FILE_EXTERNAL_INPUT_item(nid)%file_current)), &
                                                    ', varname = ', trim(FILE_EXTERNAL_INPUT_item(nid)%varname)
                      call PRC_abort
                   endif

                   FILE_EXTERNAL_INPUT_item(nid)%fid              = fid
                   FILE_EXTERNAL_INPUT_item(nid)%data_step_offset = FILE_EXTERNAL_INPUT_item(nid)%step_num
                   FILE_EXTERNAL_INPUT_item(nid)%step_num         = FILE_EXTERNAL_INPUT_item(nid)%step_num + step_nmax

                else
                   LOG_ERROR("FILE_EXTERNAL_INPUT_time_advance",*) 'Current time is out of period of external data! '
                   call PRC_abort
                endif

             endif

          endif ! periodic or not

       endif ! read from file?

    endif ! do read?


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


  function FILE_EXTERNAL_INPUT_getId( varname ) result( nid )
    character(len=*), intent(in) :: varname

    integer nid, n

    ! searching the data ID
    nid = 0
    do n = 1, FILE_EXTERNAL_INPUT_item_count
       if( varname == FILE_EXTERNAL_INPUT_item(n)%varname ) then
          nid = n
          exit
       end if
    enddo

    return
  end function FILE_EXTERNAL_INPUT_getId

end module scale_file_external_input
