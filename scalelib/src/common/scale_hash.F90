!-------------------------------------------------------------------------------
!> module CONSTANT
!!
!! @par Description
!!          Hash table module
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_hash
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_prc, only: &
     myrank => PRC_myrank,  &
     PRC_abort
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  interface hash_table
     module procedure :: table_new
  end interface hash_table
  public :: table_new

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type :: hash_entry
#if defined(__GFORTRAN__) && __GNUC__ < 7
     character(len=128), pointer :: key
#else
     character(len=:), pointer :: key
#endif
     class(*), pointer :: val
     integer :: cnt = 0
     integer :: hash
     type(hash_entry), pointer :: next => null()
  end type hash_entry

  type :: hash_entry_ptr
     type(hash_entry), pointer :: ptr => null()
  end type hash_entry_ptr

  type, public :: hash_table
     type(hash_entry_ptr), allocatable :: table(:)
     integer :: size
     integer :: len
     integer :: max_len
   contains
     procedure :: destroy => destroy
     procedure :: length => length
     procedure :: has_key => has_key
     procedure :: max_key_len => max_key_len
     procedure :: keys => keys
     procedure :: get => get
     procedure :: get_with_cnt => get_with_cnt
     procedure :: put => put
     procedure :: accumurate => accumurate
     procedure :: debug => debug
  end type hash_table

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: INIT_SIZE = 64
  integer, parameter :: HASH_P = 19
  integer, parameter :: HASH_MAX = 124901413
  integer, parameter :: SIZE_MAX = 2**26 ! 26 = floor(log_2(HASH_MAX))

contains

  !-----------------------------------------------------------------------------
  !> Constructor
  type(hash_table) function table_new() result(tbl)
    allocate( tbl%table(INIT_SIZE) )
    tbl%size = INIT_SIZE
    tbl%len = 0
    tbl%max_len = 0

    return
  end function table_new

  !-----------------------------------------------------------------------------
  !> Destructor
  subroutine destroy(self)
    class(hash_table), intent(inout) :: self
    type(hash_entry), pointer :: e, e_old
    integer :: idx

    do idx = 1, self%size
       e => self%table(idx)%ptr
       do while ( associated(e) )
          deallocate( e%key, e%val )
          e_old => e
          e => e%next
          deallocate( e_old )
       end do
    end do

    deallocate( self%table )
    self%size = 0
    self%len = 0
    self%max_len = 0

    return
  end subroutine destroy

  ! methods

  !-----------------------------------------------------------------------------
  !> Get length of elements
  integer function length(self)
    class(hash_table), intent(in) :: self

    length = self%len

  end function length

  !-----------------------------------------------------------------------------
  !> Check if key exists
  logical function has_key(self, key) result(res)
    class(hash_table), intent(in) :: self
    character(len=*), intent(in) :: key

    type(hash_entry), pointer :: e
    integer :: hash
    integer :: idx

    hash = get_hash(key)
    idx = get_index(hash, self%size)
    e => self%table(idx)%ptr
    do while ( associated(e) )
       if ( e%key == key ) then
          res = .true.
          return
       end if
       e => e%next
    end do
    res = .false.
  end function has_key

  !-----------------------------------------------------------------------------
  !> Get maximam length of keys
  integer function max_key_len(self, key)
    class(hash_table), intent(in) :: self
    character(len=*), intent(in) :: key

    max_key_len = self%max_len

  end function max_key_len

  !-----------------------------------------------------------------------------
  !> Get all keys
  subroutine keys(self, ary)
    class(hash_table), intent(in) :: self
    character(len=*), intent(out) :: ary(:)
    type(hash_entry), pointer :: e
    integer :: idx, i

    if ( size(ary) < self%len ) then
       write(*,*) "size of ary is not enough: ", self%len
       call abort()
    end if

    i = 1
    do idx = 1, self%size
       e => self%table(idx)%ptr
       do while ( associated(e) )
          ary(i) = e%key
          i = i+1
          e => e%next
       end do
    end do

  end subroutine keys

  !-----------------------------------------------------------------------------
  !> Get value
  function get(self, key) result(val)
    class(hash_table), intent(in) :: self
    character(len=*), intent(in) :: key
    class(*), pointer :: val

    type(hash_entry), pointer :: e
    integer :: hash
    integer :: idx

    hash = get_hash(key)
    idx = get_index(hash, self%size)
    e => self%table(idx)%ptr
    do while ( associated(e) )
       if ( e%key == key ) then
          val => e%val
          return
       end if
       e => e%next
    end do
    nullify(val)

  end function get

  !-----------------------------------------------------------------------------
  !> Get value and count
  subroutine get_with_cnt(self, key, val, cnt)
    class(hash_table), intent(in) :: self
    character(len=*), intent(in) :: key
    class(*), pointer, intent(out) :: val
    integer, intent(out) :: cnt

    type(hash_entry), pointer :: e
    integer :: hash
    integer :: idx

    hash = get_hash(key)
    idx = get_index(hash, self%size)
    e => self%table(idx)%ptr
    do while ( associated(e) )
       if ( e%key == key ) then
          val => e%val
          cnt = e%cnt
          return
       end if
       e => e%next
    end do

  end subroutine get_with_cnt

  !-----------------------------------------------------------------------------
  !> Put value
  subroutine put(self, key, val)
    class(hash_table), intent(inout) :: self
    character(len=*), intent(in) :: key
    class(*), intent(in) :: val

    type(hash_entry), pointer :: e
    integer :: ival
    real(RP) :: rval
    integer :: hash
    integer :: idx

    hash = get_hash(key)
    idx = get_index(hash, self%size)

    ! try to find exist entry
    e => self%table(idx)%ptr
    do while ( associated(e) )
       if ( e%key == key ) then
          deallocate( e%val )
          allocate( e%val, source=val )
          return
       end if
       e => e%next
    end do

    ! if not found
    call new_entry(self, idx, hash, key, val)

    return
  end subroutine put

  !-----------------------------------------------------------------------------
  !> Add value
  subroutine accumurate(self, key, val)
    class(hash_table), intent(inout) :: self
    character(len=*), intent(in) :: key
    real(RP), intent(in) :: val

    type(hash_entry), pointer :: e
    integer :: hash
    integer :: idx

    hash = get_hash(key)
    idx = get_index(hash, self%size)

    ! try to find exist entry
    e => self%table(idx)%ptr
    do while ( associated(e) )
       if ( e%key == key ) then
          select type ( v => e%val )
          type is ( integer )
             v = v + val
          type is ( real(RP) )
             v = v + val
          class default
             write(*,*) "type is invalid"
             call abort()
          end select
          e%cnt = e%cnt + 1
          return
       end if
       e => e%next
    end do

    call new_entry(self, idx, hash, key, val)

    return
  end subroutine accumurate


  ! private methods

  function get_hash(key) result(hash)
    character(len=*), intent(in) :: key
    integer :: hash

    integer :: i

    hash = 0
    do i = 1, len_trim(key)
       hash = mod( hash * HASH_P + ichar(key(i:i)), HASH_MAX)
    end do

  end function get_hash

  function get_index(hash, size) result(idx)
    integer, intent(in) :: hash
    integer, intent(in) :: size
    integer :: idx

    idx = iand(hash, size-1) + 1
    !write(*,*)"hash", hash, "idx", idx

  end function get_index

  subroutine new_entry(self, idx, hash, key, val)
    class(hash_table), intent(inout) :: self
    integer, intent(in) :: idx, hash
    character(len=*), intent(in) :: key
    class(*), intent(in) :: val

    type(hash_entry), pointer :: e_new
    integer :: new_size
    integer :: key_len

    ! new entry
    allocate( e_new )
    key_len = len_trim(key)
#if defined(__GFORTRAN__) && __GNUC__ < 7
    if ( key_len > 128 ) then
      LOG_ERROR("new_entry",*) 'length of the key must be <=128: ', key_len
      call PRC_abort
    endif
#endif
    allocate( e_new%key, source=key )
    allocate( e_new%val, source=val )
    e_new%cnt = 1
    e_new%hash = hash
    e_new%next => self%table(idx)%ptr
    self%table(idx)%ptr => e_new
    self%len = self%len + 1
    self%max_len = max(self%max_len, key_len)

    ! resize
    if ( self%len * 2 >= self%size ) then
       if ( self%size == SIZE_MAX ) return
       new_size = min( self%size * 2, SIZE_MAX )
       block
         type(hash_entry_ptr) :: table(new_size)
         type(hash_entry), pointer :: e, next
         integer :: i, idx_new

         do i = 1, self%size
            e => self%table(i)%ptr
            do while ( associated(e) )
               next => e%next
               idx_new = get_index(e%hash, new_size)
               e%next => table(idx_new)%ptr
               table(idx_new)%ptr => e
               e => next
            end do
         end do
         deallocate( self%table )
         allocate( self%table(new_size) )
         do i = 1, new_size
            self%table(i)%ptr => table(i)%ptr
         end do
         self%size = new_size
       end block
    end if

    return
  end subroutine new_entry

  subroutine debug(self)
    class(hash_table), intent(in) :: self
    type(hash_entry), pointer :: e
    integer :: idx, i

    write(*,'(a5,a32,a3,a8,a15,a2,a4,a10)') 'id', "key", " ", "val", "cnt", ", ", "idx", "hash"
    i = 1
    do idx = 1, self%size
       e => self%table(idx)%ptr
       do while ( associated(e) )
          select type(v => e%val)
          type is ( real(RP) )
             write(*,'(i5,a32,a3,f15.3,i8,a2,i4,i10)') i, trim(e%key), " => ", v, e%cnt, ", ", idx, e%hash
          type is ( integer )
             write(*,'(i5,a32,a3,i15,i8,a2,i4,i10)') i, trim(e%key), " => ", v, e%cnt, ", ", idx, e%hash
          class default
             write(*,'(i5,a32,a3,a15,i8,a2,i4,i10)') i, trim(e%key), " => ", "type", e%cnt, ", ", idx, e%hash
          end select
          i = i+1
          e => e%next
       end do
    end do

  end subroutine debug

end module scale_hash
