#ifndef GADG_RESERVE_WORK_X
  COMPILE_ERROR "'GADG_RESERVE_WORK_X' is not defined"
#endif
#ifndef GADG_DEALLOC_WORK_X
  COMPILE_ERROR "'GADG_DEALLOC_WORK_X' is not defined"
#endif
#ifndef GADG_WORK_X
  COMPILE_ERROR "'GADG_WORK_X' is not defined"
#endif

  !###########################################################################
  subroutine GADG_RESERVE_WORK_X(siz, stat)
    ! If the work array GADG_XWORK is not allocated or
    ! the size of GADG_XWORK is smaller than siz,
    ! GADG_XWORK is (re)allocated.
    
    integer, intent(in) :: siz
    integer, intent(out), optional :: stat
    
    !=------------------------------------------------------------------------
    if (allocated(GADG_WORK_X)) then
      if (size(GADG_WORK_X) >= siz) return
      call GADG_DEALLOC_WORK_X
    end if
    
    if (present(stat)) then
      allocate (GADG_WORK_X(siz), stat = stat)
    else
      allocate (GADG_WORK_X(siz))
    end if
    
    return
  end subroutine GADG_RESERVE_WORK_X
  
  !###########################################################################
  subroutine GADG_DEALLOC_WORK_X
    ! Deallocates the work array GADG_WORK_X.
    ! Does nothing when GADG_WORK_X is not allocated.
    
    !=------------------------------------------------------------------------
    
    if (allocated(GADG_WORK_X)) then
      deallocate (GADG_WORK_X)
    end if
    
    return
  end subroutine GADG_DEALLOC_WORK_X

#undef GADG_RESERVE_WORK_X
#undef GADG_DEALLOC_WORK_X
#undef GADG_WORK_X
