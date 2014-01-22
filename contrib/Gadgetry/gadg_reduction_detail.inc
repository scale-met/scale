# ifndef GM_SUBNAME_SEQ
  COMPILE_ERROR "GM_SUBNAME_SEQ is not defined"
# endif
# ifndef GM_SUBNAME_SPLIT
  COMPILE_ERROR "GM_SUBNAME_SPLIT is not defined"
# endif
# ifndef GM_DIM_REDUCED
  COMPILE_ERROR "GM_DIM_REDUCED is not defined"
# endif
# ifndef GM_REDUCE_OP
  COMPILE_ERROR "GM_REDUCE_OP is not defined"
# endif
# ifndef GM_INIT_R
  COMPILE_ERROR "GM_INIT_R is not defined"
# endif
# ifndef GM_INIT_V
  COMPILE_ERROR "GM_INIT_V is not defined"
# endif
# ifndef GM_START_I
  COMPILE_ERROR "GM_START_I is not defined"
# endif
# ifndef GTYPE1
  COMPILE_ERROR "GTYPE1 is not defined"
# endif
# ifndef GTYPE2
  COMPILE_ERROR "GTYPE2 is not defined"
# endif

# if GM_DIM_REDUCED == 1
#   define GM_A(I,J)  a(I,J)
#   define GM_DIM_SCAN  2
# else
#   define GM_A(I,J)  a(J,I)
#   define GM_DIM_SCAN  1
#endif

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  subroutine GM_SUBNAME_SEQ(a, rslt)
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    
    !=------------------------------------------------------------------------
    GTYPE2 :: r
    integer :: iend_reduced, iend_scan
    integer :: i, j
#ifdef GM_NEED_Q
    GTYPE1 :: q
#endif
    
    iend_reduced = ubound(a, GM_DIM_REDUCED)
    iend_scan = ubound(a, GM_DIM_SCAN)

!CDIR NOLOOPCHG
    do j = 1, iend_scan
      GM_INIT_R
      do i = GM_START_I, iend_reduced
        GM_REDUCE_OP(r, GM_A(i, j), q)
      end do
      rslt(j) = r
    end do
  end subroutine GM_SUBNAME_SEQ
  
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  subroutine GM_SUBNAME_SPLIT(a, rslt)
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    
    !=------------------------------------------------------------------------
    GTYPE2 :: v(1 : GADG_CONF_VEC_REG_LEN)
    integer :: iend_reduced, iend_scan, iend_vloop
    integer :: i, j, k
#ifdef GM_NEED_Q
    GTYPE1 :: q(1 : GADG_CONF_VEC_REG_LEN)
!CDIR VREG(q)
#endif
!CDIR VREG(v)
    
    iend_reduced = ubound(a, GM_DIM_REDUCED)
    iend_scan = ubound(a, GM_DIM_SCAN)
    
    do k = 1, iend_scan, GADG_CONF_VEC_REG_LEN
      iend_vloop = min(iend_scan - k + 1, GADG_CONF_VEC_REG_LEN)
      GM_INIT_V
      do i = GM_START_I, iend_reduced
        do j = 1, iend_vloop
          GM_REDUCE_OP(v(j), GM_A(i, j + k - 1), q(j))
        end do
      end do
      rslt(k : k + iend_vloop - 1) = v(1 : iend_vloop)
    end do
  end subroutine GM_SUBNAME_SPLIT
  
#undef GM_DIM_SCAN
#undef GM_A
