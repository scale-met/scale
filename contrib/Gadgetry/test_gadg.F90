#include "gadg_macros.inc"
module gadg_test
  use gadg_base
  use gadg_algorithm
  use gadg_reduction
  integer :: error_count=0
!  error_count = 0
 
!  call test_gadg_macros
!  call test_gadg_base
!  call test_gadg_algorithm
!  call test_gadg_reduction
!  call gadg_finalize
!  print *, "error_count = ", error_count
  
contains
  
  subroutine test_gadg_macros()
    !#= test GADG_IS_ODD(3)
    /*||*/  if (GADG_IS_ODD(3)) then
    /*||*/    print *, "  OK ... test_gadg.F90, 19"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 19"
    /*||*/    print *, 'Required: GADG_IS_ODD(3)'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test .not.GADG_IS_ODD(6)
    /*||*/  if (.not.GADG_IS_ODD(6)) then
    /*||*/    print *, "  OK ... test_gadg.F90, 27"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 27"
    /*||*/    print *, 'Required: .not.GADG_IS_ODD(6)'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test .not.GADG_IS_EVEN(3)
    /*||*/  if (.not.GADG_IS_EVEN(3)) then
    /*||*/    print *, "  OK ... test_gadg.F90, 35"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 35"
    /*||*/    print *, 'Required: .not.GADG_IS_EVEN(3)'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test GADG_IS_EVEN(6)
    /*||*/  if (GADG_IS_EVEN(6)) then
    /*||*/    print *, "  OK ... test_gadg.F90, 43"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 43"
    /*||*/    print *, 'Required: GADG_IS_EVEN(6)'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
  end subroutine test_gadg_macros
  
  
  subroutine test_gadg_base()
    call gadg_reserve_work_i32(50)
    call gadg_reserve_work_i32(20)
    !#= test [size(GADG_WORK_I32)] == 50
    /*||*/  if (size(GADG_WORK_I32) == 50) then
    /*||*/    print *, "  OK ... test_gadg.F90, 57"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 57"
    /*||*/    print *, 'size(GADG_WORK_I32) = ', size(GADG_WORK_I32)
    /*||*/    print *, 'Required: size(GADG_WORK_I32) == 50'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
  end subroutine test_gadg_base
  
  
  subroutine test_gadg_algorithm()
    real(dp) :: r
    integer, parameter :: n = 100000, nd = 200
    integer, allocatable :: key(:), freq(:), sorted(:), tag(:)
    integer, parameter :: VL = GADG_CONF_VEC_REG_LEN
    integer :: i, j, k
    allocate(key(n), freq(0 : n/nd-1), sorted(n), tag(0 : n/nd-1))
    ! gadg_vmat_height
    !#= test [gadg_vmat_height(0)] == 0
    /*||*/  if (gadg_vmat_height(0) == 0) then
    /*||*/    print *, "  OK ... test_gadg.F90, 77"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 77"
    /*||*/    print *, 'gadg_vmat_height(0) = ', gadg_vmat_height(0)
    /*||*/    print *, 'Required: gadg_vmat_height(0) == 0'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test [gadg_vmat_height(1)] == 0
    /*||*/  if (gadg_vmat_height(1) == 0) then
    /*||*/    print *, "  OK ... test_gadg.F90, 86"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 86"
    /*||*/    print *, 'gadg_vmat_height(1) = ', gadg_vmat_height(1)
    /*||*/    print *, 'Required: gadg_vmat_height(1) == 0'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test [gadg_vmat_height(VL-1)] == 0
    /*||*/  if (gadg_vmat_height(VL-1) == 0) then
    /*||*/    print *, "  OK ... test_gadg.F90, 95"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 95"
    /*||*/    print *, 'gadg_vmat_height(VL-1) = ', gadg_vmat_height(VL-1)
    /*||*/    print *, 'Required: gadg_vmat_height(VL-1) == 0'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    
    !#= test [gadg_vmat_height(VL)] == 1
    /*||*/  if (gadg_vmat_height(VL) == 1) then
    /*||*/    print *, "  OK ... test_gadg.F90, 105"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 105"
    /*||*/    print *, 'gadg_vmat_height(VL) = ', gadg_vmat_height(VL)
    /*||*/    print *, 'Required: gadg_vmat_height(VL) == 1'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test [gadg_vmat_height(VL+1)] == 1
    /*||*/  if (gadg_vmat_height(VL+1) == 1) then
    /*||*/    print *, "  OK ... test_gadg.F90, 114"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 114"
    /*||*/    print *, 'gadg_vmat_height(VL+1) = ', gadg_vmat_height(VL+1)
    /*||*/    print *, 'Required: gadg_vmat_height(VL+1) == 1'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test [gadg_vmat_height(VL*2-1)] == 1
    /*||*/  if (gadg_vmat_height(VL*2-1) == 1) then
    /*||*/    print *, "  OK ... test_gadg.F90, 123"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 123"
    /*||*/    print *, 'gadg_vmat_height(VL*2-1) = ', gadg_vmat_height(VL*2-1)
    /*||*/    print *, 'Required: gadg_vmat_height(VL*2-1) == 1'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    
    !#= test [gadg_vmat_height(VL*2)] == 1
    /*||*/  if (gadg_vmat_height(VL*2) == 1) then
    /*||*/    print *, "  OK ... test_gadg.F90, 133"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 133"
    /*||*/    print *, 'gadg_vmat_height(VL*2) = ', gadg_vmat_height(VL*2)
    /*||*/    print *, 'Required: gadg_vmat_height(VL*2) == 1'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test [gadg_vmat_height(VL*2+1)] == 1
    /*||*/  if (gadg_vmat_height(VL*2+1) == 1) then
    /*||*/    print *, "  OK ... test_gadg.F90, 142"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 142"
    /*||*/    print *, 'gadg_vmat_height(VL*2+1) = ', gadg_vmat_height(VL*2+1)
    /*||*/    print *, 'Required: gadg_vmat_height(VL*2+1) == 1'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test [gadg_vmat_height(VL*3-1)] == 1
    /*||*/  if (gadg_vmat_height(VL*3-1) == 1) then
    /*||*/    print *, "  OK ... test_gadg.F90, 151"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 151"
    /*||*/    print *, 'gadg_vmat_height(VL*3-1) = ', gadg_vmat_height(VL*3-1)
    /*||*/    print *, 'Required: gadg_vmat_height(VL*3-1) == 1'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    
    !#= test [gadg_vmat_height(VL*3)] == 3
    /*||*/  if (gadg_vmat_height(VL*3) == 3) then
    /*||*/    print *, "  OK ... test_gadg.F90, 161"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 161"
    /*||*/    print *, 'gadg_vmat_height(VL*3) = ', gadg_vmat_height(VL*3)
    /*||*/    print *, 'Required: gadg_vmat_height(VL*3) == 3'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test [gadg_vmat_height(VL*3+1)] == 3
    /*||*/  if (gadg_vmat_height(VL*3+1) == 3) then
    /*||*/    print *, "  OK ... test_gadg.F90, 170"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 170"
    /*||*/    print *, 'gadg_vmat_height(VL*3+1) = ', gadg_vmat_height(VL*3+1)
    /*||*/    print *, 'Required: gadg_vmat_height(VL*3+1) == 3'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    !#= test [gadg_vmat_height(VL*4-1)] == 3
    /*||*/  if (gadg_vmat_height(VL*4-1) == 3) then
    /*||*/    print *, "  OK ... test_gadg.F90, 179"
    /*||*/  else
    /*||*/    print *, "Test fails! ... test_gadg.F90, 179"
    /*||*/    print *, 'gadg_vmat_height(VL*4-1) = ', gadg_vmat_height(VL*4-1)
    /*||*/    print *, 'Required: gadg_vmat_height(VL*4-1) == 3'
    /*||*/    error_count = error_count + 1
    /*||*/  end if
    
    ! gadg_accumulate
    call random_seed()
    do i = 1, n
      call random_number(r)
      sorted(i) = int(r * 10000)
    end do
    call gadg_accumulate(sorted, key)
    k = 0
    do i = 1, n
      k = k + sorted(i)
      !#= qtest [k] == [key(i)]
      /*||*/  if (.not.(k == key(i))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 199"
      /*||*/    print *, 'k = ', k
      /*||*/    print *, 'key(i) = ', key(i)
      /*||*/    print *, 'Required: k == key(i)'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_accumulate(sorted)
    do i = 1, n
      !#= qtest [sorted(i)] == [key([i])]
      /*||*/  if (.not.(sorted(i) == key(i))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 210"
      /*||*/    print *, 'sorted(i) = ', sorted(i)
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'key(i) = ', key(i)
      /*||*/    print *, 'Required: sorted(i) == key(i)'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    
    ! gadg_count_sort
    call random_seed()
    do i = 0, n / nd - 1
      do j = 0, nd - 1
        key(i + 1 + j * (n / nd)) = i
      end do
    end do
    do i = 1, n - 1
      call random_number(r)
      j = int(r * (n - i)) + i + 1
      k = key(i)
      key(i) = key(j)
      key(j) = k
    end do
    
    call gadg_count_sort(key, 0, n/nd-1, freq, tag, sorted, local_work = .true.)
    do i = 0, n / nd - 1
      !#= qtest [freq([i])] == nd
      /*||*/  if (.not.(freq(i) == nd)) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 237"
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'freq(i) = ', freq(i)
      /*||*/    print *, 'Required: freq(i) == nd'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
      !#= qtest [tag([i])] == [i * nd + 1]
      /*||*/  if (.not.(tag(i) == i * nd + 1)) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 245"
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'tag(i) = ', tag(i)
      /*||*/    print *, 'i * nd + 1 = ', i * nd + 1
      /*||*/    print *, 'Required: tag(i) == i * nd + 1'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    do i = 1, n
      !#= qtest [key([sorted([i])])] == (i - 1) / nd
      /*||*/  if (.not.(key(sorted(i)) == (i - 1) / nd)) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 256"
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'sorted(i) = ', sorted(i)
      /*||*/    print *, 'key(sorted(i)) = ', key(sorted(i))
      /*||*/    print *, 'Required: key(sorted(i)) == (i - 1) / nd'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
      key(sorted(i)) = -1
    end do
  end subroutine test_gadg_algorithm
  
  
  subroutine test_gadg_reduction()
    real(dp) :: r
    integer :: i, j
    integer, allocatable :: a(:,:), rslt(:)
    integer, parameter :: m = 70, n = 10001
    
    allocate(a(m,n), rslt(n))
    do i = 1, m
      do j = 1, n
        call random_number(r)
        a(i,j) = int(r * 100)
      end do
    end do
    call gadg_reduction_minval(a, rslt, 1)
    do j = 1, n
      !#= qtest [rslt([j])] == [minval(a(:,j))]
      /*||*/  if (.not.(rslt(j) == minval(a(:,j)))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 285"
      /*||*/    print *, 'j = ', j
      /*||*/    print *, 'rslt(j) = ', rslt(j)
      /*||*/    print *, 'minval(a(:,j)) = ', minval(a(:,j))
      /*||*/    print *, 'Required: rslt(j) == minval(a(:,j))'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_reduction_maxval(a, rslt, 1)
    do j = 1, n
      !#= qtest [rslt([j])] == [maxval(a(:,j))]
      /*||*/  if (.not.(rslt(j) == maxval(a(:,j)))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 297"
      /*||*/    print *, 'j = ', j
      /*||*/    print *, 'rslt(j) = ', rslt(j)
      /*||*/    print *, 'maxval(a(:,j)) = ', maxval(a(:,j))
      /*||*/    print *, 'Required: rslt(j) == maxval(a(:,j))'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_reduction_minloc(a, rslt, 1)
    do j = 1, n
      !#= qtest [rslt([j])] == [minloc(a(:,j),dim=1)]
      /*||*/  if (.not.(rslt(j) == minloc(a(:,j),dim=1))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 309"
      /*||*/    print *, 'j = ', j
      /*||*/    print *, 'rslt(j) = ', rslt(j)
      /*||*/    print *, 'minloc(a(:,j),dim=1) = ', minloc(a(:,j),dim=1)
      /*||*/    print *, 'Required: rslt(j) == minloc(a(:,j),dim=1)'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_reduction_maxloc(a, rslt, 1)
    do j = 1, n
      !#= qtest [rslt([j])] == [maxloc(a(:,j),dim=1)]
      /*||*/  if (.not.(rslt(j) == maxloc(a(:,j),dim=1))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 321"
      /*||*/    print *, 'j = ', j
      /*||*/    print *, 'rslt(j) = ', rslt(j)
      /*||*/    print *, 'maxloc(a(:,j),dim=1) = ', maxloc(a(:,j),dim=1)
      /*||*/    print *, 'Required: rslt(j) == maxloc(a(:,j),dim=1)'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_reduction_sum(a, rslt, 1)
    do j = 1, n
      !#= qtest [rslt([j])] == [sum(a(:,j))]
      /*||*/  if (.not.(rslt(j) == sum(a(:,j)))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 333"
      /*||*/    print *, 'j = ', j
      /*||*/    print *, 'rslt(j) = ', rslt(j)
      /*||*/    print *, 'sum(a(:,j)) = ', sum(a(:,j))
      /*||*/    print *, 'Required: rslt(j) == sum(a(:,j))'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    
    call gadg_reduction_minval(a, rslt, 2)
    do i = 1, m
      !#= qtest [rslt([i])] == [minval(a(i,:))]
      /*||*/  if (.not.(rslt(i) == minval(a(i,:)))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 346"
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'rslt(i) = ', rslt(i)
      /*||*/    print *, 'minval(a(i,:)) = ', minval(a(i,:))
      /*||*/    print *, 'Required: rslt(i) == minval(a(i,:))'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_reduction_maxval(a, rslt, 2)
    do i = 1, m
      !#= qtest [rslt([i])] == [maxval(a(i,:))]
      /*||*/  if (.not.(rslt(i) == maxval(a(i,:)))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 358"
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'rslt(i) = ', rslt(i)
      /*||*/    print *, 'maxval(a(i,:)) = ', maxval(a(i,:))
      /*||*/    print *, 'Required: rslt(i) == maxval(a(i,:))'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_reduction_minloc(a, rslt, 2)
    do i = 1, m
      !#= qtest [rslt([i])] == [minloc(a(i,:),dim=1)]
      /*||*/  if (.not.(rslt(i) == minloc(a(i,:),dim=1))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 370"
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'rslt(i) = ', rslt(i)
      /*||*/    print *, 'minloc(a(i,:),dim=1) = ', minloc(a(i,:),dim=1)
      /*||*/    print *, 'Required: rslt(i) == minloc(a(i,:),dim=1)'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_reduction_maxloc(a, rslt, 2)
    do i = 1, m
      !#= qtest [rslt([i])] == [maxloc(a(i,:),dim=1)]
      /*||*/  if (.not.(rslt(i) == maxloc(a(i,:),dim=1))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 382"
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'rslt(i) = ', rslt(i)
      /*||*/    print *, 'maxloc(a(i,:),dim=1) = ', maxloc(a(i,:),dim=1)
      /*||*/    print *, 'Required: rslt(i) == maxloc(a(i,:),dim=1)'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    call gadg_reduction_sum(a, rslt, 2)
    do i = 1, m
      !#= qtest [rslt([i])] == [sum(a(i,:))]
      /*||*/  if (.not.(rslt(i) == sum(a(i,:)))) then
      /*||*/    print *, "Test fails! ... test_gadg.F90, 394"
      /*||*/    print *, 'i = ', i
      /*||*/    print *, 'rslt(i) = ', rslt(i)
      /*||*/    print *, 'sum(a(i,:)) = ', sum(a(i,:))
      /*||*/    print *, 'Required: rslt(i) == sum(a(i,:))'
      /*||*/    error_count = error_count + 1
      /*||*/  end if
    end do
    
  end subroutine test_gadg_reduction
  
end module gadg_test
