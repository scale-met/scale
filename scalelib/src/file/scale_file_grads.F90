!-------------------------------------------------------------------------------
!> module file_grads
!!
!! @par Description
!!          read data from GrADS file
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_file_grads
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: FILE_GrADS_open
  public :: FILE_GrADS_varid
  public :: FILE_GrADS_isOneD
  public :: FILE_GrADS_get_shape
  public :: FILE_GrADS_read
  public :: FILE_GrADS_close

  interface FILE_GrADS_get_shape
     module procedure FILE_GrADS_get_shape_name
     module procedure FILE_GrADS_get_shape_id
  end interface FILE_GrADS_get_shape

  interface FILE_GrADS_read
     module procedure FILE_GrADS_read_1D_name
     module procedure FILE_GrADS_read_1D_id
     module procedure FILE_GrADS_read_2D_name
     module procedure FILE_GrADS_read_2D_id
     module procedure FILE_GrADS_read_3D_name
     module procedure FILE_GrADS_read_3D_id
  end interface FILE_GrADS_read
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
  integer, parameter :: nmls_max = 10
  integer, parameter :: vars_max  = 100
  integer, parameter :: lvars_max = 1000
  type t_var
     character(len=H_SHORT) :: name
     character(len=H_LONG)  :: fname
     character(len=H_SHORT) :: dtype
     real(RP)               :: swpoint
     real(RP)               :: dd
     integer                :: lnum
     real(RP), allocatable  :: lvars(:)
     integer                :: startrec
     integer                :: totalrec
     real(SP)               :: missval
     integer                :: nx
     integer                :: ny
     integer                :: nz
     logical                :: yrev
     integer                :: endian ! 0: little, 1: big
  end type t_var
  type t_nml
     character(len=H_LONG)    :: fname
     integer                  :: nx, ny, nz
     type(t_var), allocatable :: vars(:)
     integer                  :: nvars
  end type t_nml
  type t_file
     character(len=H_LONG)  :: fname
     character(len=H_SHORT) :: postfix
     integer                :: fid
  end type t_file
  type(t_nml), save  :: nmls(nmls_max)
  integer      :: nnmls = 0
  type(t_file) :: files(vars_max)
  integer      :: nfiles = 0

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Open
  subroutine FILE_GrADS_open( &
       file_name, &
       file_id    )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    character(len=*), intent(in)  :: file_name
    integer,          intent(out) :: file_id

    character(len=H_SHORT) :: name             ! up to 16 characters
    character(len=H_SHORT) :: dtype            ! 'linear','levels','map'
    character(len=H_LONG)  :: fname            ! head of file name
    real(RP)               :: swpoint          ! start point (south-west point) for "linear"
    real(RP)               :: dd               ! dlon,dlat for "linear"
    integer                :: lnum             ! number of data
    real(RP)               :: lvars(lvars_max) ! values for "levels"
    integer                :: startrec         ! record position
    integer                :: totalrec         ! total record number per one time
    real(SP)               :: missval          ! missing value
    integer                :: nx               ! optional
    integer                :: ny               ! optional
    integer                :: nz               ! optional
    character(len=H_SHORT) :: fendian          ! option for "map"
    logical                :: yrev             ! option for "map", if yrev=.true., order of data is NW to SE.

    namelist /GrADS_DIMS/ &
       nx, &
       ny, &
       nz

    namelist /GrADS_ITEM/ &
       name,     & ! necessary
       dtype,    & ! necessary
       fname,    & ! necessary except for linear data
       swpoint,  & ! for linear data
       dd,       & ! for linear data
       lnum,     & ! for levels data
       lvars,    & ! for levels data
       startrec, & ! for map data
       totalrec, & ! for map data
       missval,  & ! option
       nx,       & ! option
       ny,       & ! option
       nz,       & ! option
       yrev        ! option
!       fendian     ! option

    integer :: fid
    integer :: nvars
    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("FILE_GrADS_open",*) 'open namelist file :', trim(file_name)

    ! check exist
    do n = 1, nnmls
       if ( nmls(n)%fname == file_name ) then
          ! alread read
          file_id = n
          return
       end if
    end do


    fid = IO_get_available_fid()
    !--- open namelist
    open( fid, &
          file   = trim(file_name), &
          form   = 'formatted',    &
          status = 'old',          &
          action = 'read',         &
          iostat = ierr            )
    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_GrADS_open",*) 'Input file is not found! ', trim(file_name)
       call PRC_abort
    endif

    call check_oldnamelist( fid )

    !--- read namelist dims
    read(fid,nml=GrADS_DIMS,iostat=ierr)
    if( ierr /= 0 ) then !--- missing or fatal error
       LOG_ERROR("FILE_GrADS_open",*) 'Not appropriate names in GrADS_DIMS in ', trim(file_name),'. Check!'
       call PRC_abort
    endif
    LOG_NML(GrADS_DIMS)


    nnmls = nnmls + 1
    file_id = nnmls

    nmls(file_id)%fname = file_name
    nmls(file_id)%nx    = nx
    nmls(file_id)%ny    = ny
    nmls(file_id)%nz    = nz



    !--- count the number of variables
    rewind(fid)
    nvars = 0
    do n = 1, vars_max
       read(fid, nml=GrADS_ITEM, iostat=ierr)
       if( ierr > 0 )then
          LOG_ERROR("FILE_GrADS_open",*) 'Not appropriate names in GrADS_ITEM in ', trim(file_name),'. Check!'
          call PRC_abort
       else if( ierr < 0 )then
          exit
       endif
       nvars = nvars + 1
    enddo

    if ( nvars > vars_max ) then
       LOG_ERROR("FILE_GRADS_open",*) 'The number of grads vars exceeds the limit! ', &
            nvars, ' > ', vars_max
       call PRC_abort
    endif

    nmls(file_id)%nvars = nvars
    allocate( nmls(file_id)%vars(nvars) )

    !--- read information of the variables
    rewind(fid)
    do n = 1, nvars

       ! set default
       name     = ''
       dtype    = ''
       fname    = ''
       swpoint  = UNDEF
       dd       = UNDEF
       lnum     = -1
       lvars(:) = UNDEF
       startrec = -1
       totalrec = -1
       nx       = nmls(file_id)%nx
       ny       = nmls(file_id)%ny
       nz       = nmls(file_id)%nz
       yrev     = .false.
       fendian  = 'big'
       missval  = UNDEF

       ! read namelist
       read(fid, nml=GrADS_ITEM, iostat=ierr)
       if( ierr /= 0 ) exit

       nmls(file_id)%vars(n)%name    = name
       nmls(file_id)%vars(n)%fname   = fname
       nmls(file_id)%vars(n)%dtype   = dtype
       nmls(file_id)%vars(n)%swpoint = swpoint
       nmls(file_id)%vars(n)%dd      = dd
       nmls(file_id)%vars(n)%lnum    = lnum
       if ( lnum > 0 ) then
          allocate( nmls(file_id)%vars(n)%lvars(lnum) )
          nmls(file_id)%vars(n)%lvars(:) = lvars(1:lnum)
       end if
       nmls(file_id)%vars(n)%startrec = startrec
       nmls(file_id)%vars(n)%totalrec = totalrec
       nmls(file_id)%vars(n)%missval  = missval
       nmls(file_id)%vars(n)%nx       = nx
       nmls(file_id)%vars(n)%ny       = ny
       nmls(file_id)%vars(n)%nz       = nz
       nmls(file_id)%vars(n)%yrev     = yrev
       if ( fendian == "big" ) then
          nmls(file_id)%vars(n)%endian = 1
       else
          nmls(file_id)%vars(n)%endian = 0
       end if

    end do

    close( fid )

  end subroutine FILE_GrADS_open

  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_varid( &
       file_id,  &
       var_name, &
       var_id    )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(in)  :: file_id
    character(len=*), intent(in)  :: var_name
    integer,          intent(out) :: var_id

    integer :: n

    if ( file_id < 0 ) then
       LOG_ERROR("FILE_GrADS_varid",*) 'file_id is invalid: ', file_id
       call PRC_abort
    end if

    var_id = -1
    do n = 1, nmls(file_id)%nvars
       if ( nmls(file_id)%vars(n)%name == var_name ) then
          var_id = n
          return
       end if
    end do

    return
  end subroutine FILE_GrADS_varid

  function FILE_GrADS_isOneD( &
       file_id,  &
       var_id ) &
       result(ret)
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in)  :: file_id
    integer, intent(in)  :: var_id
    logical :: ret

    if ( file_id < 0 ) then
       LOG_ERROR("FILE_GrADS_isOneD",*) 'file_id is invalid: ', file_id
       call PRC_abort
    end if
    if ( var_id < 0 ) then
       LOG_ERROR("FILE_GrADS_isOneD",*) 'var_id is invalid: ', var_id
       call PRC_abort
    end if

    select case( nmls(file_id)%vars(var_id)%dtype )
    case ('linear', 'levels')
       ret = .true.
    case default
       ret = .false.
    end select

    return
  end function FILE_GrADS_isOneD

  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_get_shape_name( &
       file_id,  &
       var_name, &
       shape     )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(in)  :: file_id
    character(len=*), intent(in)  :: var_name
    integer,          intent(out) :: shape(:)
    
    integer :: var_id

    call FILE_GrADS_varid( file_id, var_name, & ! (in)
                           var_id             ) ! (out)
    
    if ( var_id < 0 ) then
       LOG_ERROR("FILE_GrADS_get_shape",*) 'variable "', trim(var_name), ' is not founed in file "', trim(nmls(file_id)%fname), '"'
       call PRC_abort
    end if

    call FILE_GrADS_get_shape_id( file_id, var_id, & ! (in)
                                  shape(:)         ) ! (out)

    return
  end subroutine FILE_GrADS_get_shape_name
  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_get_shape_id( &
       file_id, &
       var_id,  &
       shape    )
    implicit none
    integer, intent(in)  :: file_id
    integer, intent(in)  :: var_id
    integer, intent(out) :: shape(:)
    
    intrinsic :: size

    if ( FILE_GrADS_isOneD( file_id, var_id ) ) then
       if ( nmls(file_id)%vars(var_id)%dtype == "levels" ) then
          shape(1) = nmls(file_id)%vars(var_id)%lnum
       else
          shape(1) = -1
       end if
    else if ( size(shape) == 2 ) then
       shape(1) = nmls(file_id)%vars(var_id)%nx
       shape(2) = nmls(file_id)%vars(var_id)%ny
    else
       shape(1) = nmls(file_id)%vars(var_id)%nz
       shape(2) = nmls(file_id)%vars(var_id)%nx
       shape(3) = nmls(file_id)%vars(var_id)%ny
    end if

    return
  end subroutine FILE_GrADS_get_shape_id

  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_read_1D_name( &
       file_id,  &
       var_name, &
       var,      &
       postfix   )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(in)  :: file_id
    character(len=*), intent(in)  :: var_name
    real(RP),         intent(out) :: var(:)
    character(len=*), intent(in), optional :: postfix

    integer :: var_id
    !---------------------------------------------------------------------------

    call FILE_GrADS_varid( file_id, var_name, & ! (in)
                           var_id             ) ! (out)

    if ( var_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_1D_name",*) 'variable "', trim(var_name), ' is not founed in file "', trim(nmls(file_id)%fname), '"'
       call PRC_abort
    end if


    call FILE_GrADS_read_1D_id( file_id, var_id,   & ! (in)
                                var(:),            & ! (out)
                                postfix = postfix  ) ! (in)

    return
  end subroutine FILE_GrADS_read_1D_name
  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_read_1D_id( &
       file_id, &
       var_id,  &
       var,     &
       postfix  )
    implicit none
    integer,          intent(in)  :: file_id
    integer,          intent(in)  :: var_id
    real(RP),         intent(out) :: var(:)
    character(len=*), intent(in), optional :: postfix

    logical :: exist
    integer :: vid

    intrinsic :: size
    !---------------------------------------------------------------------------

    if ( file_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_1D_vid",*) 'file_id is invalid: ', file_id
    end if
    if ( var_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_1D_vid",*) 'var_id is invalid: ', var_id
    end if

    call FILE_GrADS_read_data( nmls(file_id)%vars(var_id), & ! (in)
                               1, size(var),               & ! (in)
                               var(:),                     & ! (out)
                               postfix = postfix           ) ! (in)

    return
  end subroutine FILE_GrADS_read_1D_id

  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_read_2D_name( &
       file_id,  &
       var_name, &
       var,      &
       it,       &
       postfix   )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(in)  :: file_id
    character(len=*), intent(in)  :: var_name
    real(RP),         intent(out) :: var(:,:)
    integer,          intent(in), optional :: it
    character(len=*), intent(in), optional :: postfix

    integer :: var_id
    !---------------------------------------------------------------------------

    call FILE_GrADS_varid( file_id, var_name, & ! (in)
                           var_id             ) ! (out)

    if ( var_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_2D_name",*) 'variable "', trim(var_name), ' is not founed in file "', trim(nmls(file_id)%fname), '"'
       call PRC_abort
    end if

    call FILE_GrADS_read_2D_id( file_id, var_id,   & ! (in)
                                var(:,:),          & ! (out)
                                it = it,           & ! (in)
                                postfix = postfix  ) ! (in)

    return
  end subroutine FILE_GrADS_read_2D_name
  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_read_2D_id( &
       file_id, &
       var_id,  &
       var,     &
       it,      &
       postfix  )
    implicit none
    integer,          intent(in)  :: file_id
    integer,          intent(in)  :: var_id
    real(RP),         intent(out) :: var(:,:)
    integer,          intent(in), optional :: it
    character(len=*), intent(in), optional :: postfix

    integer :: vid
    intrinsic :: size
    !---------------------------------------------------------------------------

    if ( file_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_2D_vid",*) 'file_id is invalid: ', file_id
    end if
    if ( var_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_2D_vid",*) 'var_id is invalid: ', var_id
    end if

    call FILE_GrADS_read_data( nmls(file_id)%vars(var_id), & ! (in)
                               2, size(var),               & ! (in)
                               var(:,:),                   & ! (out)
                               it = it, postfix = postfix  ) ! (in)

    return
  end subroutine FILE_GrADS_read_2D_id

  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_read_3D_name( &
       file_id,  &
       var_name, &
       var,      &
       it,       &
       postfix   )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(in)  :: file_id
    character(len=*), intent(in)  :: var_name
    real(RP),         intent(out) :: var(:,:,:)
    integer,          intent(in), optional :: it
    character(len=*), intent(in), optional :: postfix

    integer :: var_id
    !---------------------------------------------------------------------------

    call FILE_GrADS_varid( file_id, var_name, & ! (in)
                           var_id             ) ! (out)

    if ( var_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_3D_name",*) 'variable "', trim(var_name), ' is not founed in file "', trim(nmls(file_id)%fname), '"'
       call PRC_abort
    end if

    call FILE_GrADS_read_3D_id( file_id, var_id,   & ! (in)
                                var(:,:,:),        & ! (out)
                                it = it,           & ! (in)
                                postfix = postfix  ) ! (in)

    return
  end subroutine FILE_GrADS_read_3D_name
  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_read_3D_id( &
       file_id, &
       var_id,  &
       var,     &
       it,      &
       postfix  )
    implicit none
    integer,          intent(in)  :: file_id
    integer,          intent(in)  :: var_id
    real(RP),         intent(out) :: var(:,:,:)
    integer,          intent(in), optional :: it
    character(len=*), intent(in), optional :: postfix

    integer :: vid
    intrinsic :: size
    !---------------------------------------------------------------------------

    if ( file_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_3D_vid",*) 'file_id is invalid: ', file_id
    end if
    if ( var_id < 0 ) then
       LOG_ERROR("FILE_GrADS_read_3D_vid",*) 'var_id is invalid: ', var_id
    end if

    call FILE_GrADS_read_data( nmls(file_id)%vars(var_id), & ! (in)
                               3, size(var),               & ! (in)
                               var(:,:,:),                 & ! (out)
                               it = it, postfix = postfix  ) ! (in)

    return
  end subroutine FILE_GrADS_read_3D_id

  !-----------------------------------------------------------------------------
  subroutine FILE_GrADS_close( &
       file_id )
    implicit none
    integer, intent(in) :: file_id

    integer :: n, m

    if ( file_id < 0 ) return

    do n = 1, nmls(file_id)%nvars
       do m = 1, nfiles
          if ( files(m)%fname == nmls(file_id)%vars(n)%fname ) then
             if ( files(m)%fid > 0 ) then
                close( files(m)%fid )
                files(m)%fid = -1
                files(m)%postfix = ""
             end if
             exit
          end if
       end do
       if ( nmls(file_id)%vars(n)%lnum > 0 ) deallocate( nmls(file_id)%vars(n)%lvars )
       nmls(file_id)%vars(n)%lnum = -1
    end do
    deallocate( nmls(file_id)%vars )
    nmls(file_id)%fname = ""
    nmls(file_id)%nvars = 0

    return
  end subroutine FILE_GrADS_close

  !-----------------------------------------------------------------------------
  ! private
  !-----------------------------------------------------------------------------

  subroutine FILE_GrADS_read_data( &
       var_info, &
       ndims, n, &
       var,      &
       it,       &
       postfix   )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       EPS   => CONST_EPS
    implicit none
    type(t_var), intent(in)  :: var_info
    integer,     intent(in)  :: ndims
    integer,     intent(in)  :: n
    real(RP),    intent(out) :: var(n)
    integer,          intent(in), optional :: it
    character(len=*), intent(in), optional :: postfix

    integer               :: fid
    character(len=H_LONG) :: gfile
    real(SP)              :: buf(var_info%nx,var_info%ny)

    integer                :: it_
    character(len=H_SHORT) :: postfix_

    integer :: nxy, nz
    integer :: irecl, ierr
    integer :: i, j, k

    select case( var_info%dtype )
    case("linear")

       if ( ndims > 1 ) then
          LOG_ERROR("FILE_GrADS_read_data",*) '"linear" is invalid for dtype of 2D or 2D var!'
       end if

       if( var_info%swpoint == UNDEF .or.var_info%dd == UNDEF ) then
          LOG_ERROR("FILE_GrADS_read_data",*) '"swpoint" and "dd" are required for linear data! ', var_info%swpoint
          call PRC_abort
       endif

       do i = 1, n
          var(i) = var_info%swpoint + (i-1) * var_info%dd
       end do

    case("levels")

       if ( ndims > 1 ) then
          LOG_ERROR("FILE_GrADS_read_data",*) '"levels" is invalid for dtype of 2D or 2D var!'
       end if

       if ( var_info%lnum < 0 )then
          LOG_ERROR("FILE_GrADS_read_data",*) '"lnum" is required for levels data! '
          call PRC_abort
       endif
       if ( var_info%lnum .ne. n ) then
          LOG_ERROR("FILE_GrADS_read_data",*) '"lnum" and size of var are not the same', var_info%lnum, n
          call PRC_abort
       end if

       do k = 1, var_info%lnum
          var(k) = var_info%lvars(k)
          if( var(k) == UNDEF ) then
             LOG_ERROR("FILE_GrADS_read_data",*) '"lvars" must be specified for levels data! '
             call PRC_abort
          endif
       end do

    case("map")

       if ( present(postfix) ) then
          postfix_ = postfix
       else
          postfix_ = ""
       end if
       if ( present(it) ) then
          it_ = it
       else
          it_ = 1
       end if

       if ( ndims == 1 ) then
          LOG_ERROR("FILE_GrADS_read_data",*) '"map" is invalid for dtype of 1D var!'
       end if

       if( var_info%startrec < 0 .or. var_info%totalrec < 0 )then
          LOG_ERROR("FILE_GrADS_read_data",*) '"startrec" and "totalrec" are required for map data! ', var_info%startrec, var_info%totalrec
          call PRC_abort
       endif
       if( var_info%fname == "" )then
          LOG_ERROR("FILE_GrADS_read_data",*) '"fname" is required for map data!'
          call PRC_abort
       endif

       ! get file_id
       fid = -1
       do i = 1, nfiles
          if ( files(i)%fname == var_info%fname ) then
             fid = i
             exit
          end if
       end do
       if ( fid < 0 ) then
          nfiles = nfiles + 1
          fid = nfiles
          files(fid)%fname   = var_info%fname
          files(fid)%postfix = ""
          files(fid)%fid     = -1
       end if

       gfile = trim(var_info%fname)//trim(postfix_)//'.grd'

       if ( files(fid)%postfix == postfix_ .and. files(fid)%fid > 0 ) then
          fid = files(fid)%fid
       else
          if ( files(fid)%fid > 0 ) close( files(fid)%fid )
          files(fid)%fid = IO_get_available_fid()
          files(fid)%postfix = postfix_
          fid = files(fid)%fid
          irecl = var_info%nx * var_info%ny * 4
          open( fid, &
                file   = gfile, &
                form   = 'unformatted', &
                access = 'direct', &
                recl   = irecl, &
                status = 'old', &
                iostat = ierr )
          if ( ierr /= 0 ) then
             LOG_ERROR("FILE_GrADS_read_data",*) 'Failed to open the grads data file! ', trim(gfile)
             call PRC_abort
          end if
       end if

       if ( ndims == 2 ) then
          nz = 1
       else
          nz = var_info%nz
       end if
       nxy = var_info%nx * var_info%ny

       if ( n .ne. nxy * nz ) then
          LOG_ERROR("FILE_GrADS_read_data",*) 'size of var is not consitent with namelist info! ', n, var_info%nx, var_info%ny, var_info%nz
          call PRC_abort
       end if

       do k = 1, nz
          irecl = var_info%totalrec * (it_-1) + var_info%startrec + k - 1
          read(fid, rec=irecl, iostat=ierr) buf(:,:)
          if ( ierr /= 0 ) then
             LOG_ERROR("FILE_GrADS_read_data",*) 'Failed to read data! ', trim(var_info%name), ', k=',k,', it=',it_, ' in ', trim(gfile)
             LOG_ERROR_CONT(*) 'irec=', irecl
             call PRC_abort
          end if

          if ( var_info%yrev ) then
             !$omp parallel do collapse(2)
             do j = 1, var_info%ny
             do i = 1, var_info%nx
                var(k+(i-1)*nz+(j-1)*var_info%nx*nz) = buf(i,var_info%ny-j+1)
             end do
             end do
          else
             !$omp parallel do collapse(2)
             do j = 1, var_info%ny
             do i = 1, var_info%nx
                var(k+(i-1)*nz+(j-1)*var_info%nx*nz) = buf(i,j)
             end do
             end do
          end if

       end do

       if ( var_info%missval .ne. UNDEF ) then
          !$omp parallel do
          do i = 1, nz * nxy
             if ( abs( var(i) - var_info%missval ) < EPS ) var(i) = UNDEF
          end do
       end if


    end select

    return
  end subroutine FILE_GrADS_read_data

  subroutine check_oldnamelist( fid )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: fid

    integer :: ierr
    logical :: dummy

    namelist /nml_grads_grid/ dummy
    namelist /grdvar/         dummy

    read(fid, nml=nml_grads_grid, iostat=ierr)
    if( ierr > 0 )then
       LOG_ERROR("check_oldnamelist",*) 'The old namelist "nml_grads_grid" is found.'
       LOG_ERROR_CONT(*) 'Use "GrADS_DIMS" instead.'
       call PRC_abort
    endif
    rewind(fid)

    read(fid, nml=grdvar, iostat=ierr)
    if( ierr > 0 )then
       LOG_ERROR("check_oldnamelist",*) 'The old namelist "grdvar" is found.'
       LOG_ERROR_CONT(*) 'Use "GrADS_ITEM" instead.'
       call PRC_abort
    endif
    rewind(fid)

    return
  end subroutine check_oldnamelist

end module scale_file_grads
