!-------------------------------------------------------------------------------
!> module scale_atmos_dyn_tinteg_rkutil
!!
!! @par Description
!!      Utility for some runge-kutta schemes
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_rkutil
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
#if defined DEBUG || defined QUICKDEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public :: RKUtil
    integer, allocatable :: comm_ind(:)
    real(RP), allocatable :: work0(:,:,:,:)
    real(RP), allocatable :: work(:,:,:,:,:)
    real(RP), allocatable :: buf(:,:,:,:)
    integer :: var_num
    integer :: register_num
  end type

  public :: ATMOS_DYN_Tinteg_RKUtil_setup
  public :: ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc
  public :: ATMOS_DYN_Tinteg_RKUtil_rkwork_dealloc
  public :: ATMOS_DYN_Tinteg_RKUtil_comm
  public :: ATMOS_DYN_Tinteg_RKUtil_comm_wait

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !----------------------------------------------------------------------------

  !- coeffecients for 7 stage RK with 6th order --------------

  real(RP), parameter, public :: RKCoef_a_7s6o_Butcher1964(7,7) = reshape( &
  (/ 0.0_RP, 1.0_RP/3.0_RP,        0.0_RP,  1.0_RP/12.0_RP, -1.0_RP/16.0_RP,         0.0_RP,   9.0_RP/44.0_RP,    &
     0.0_RP,        0.0_RP, 2.0_RP/3.0_RP,   1.0_RP/3.0_RP,   9.0_RP/8.0_RP,  9.0_RP/8.0_RP,  -9.0_RP/11.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP, -1.0_RP/12.0_RP, -3.0_RP/16.0_RP,   -3_RP/8.0_RP,  63.0_RP/44.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,  -3.0_RP/8.0_RP, -3.0_RP/4.0_RP,  18.0_RP/11.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,  1.0_RP/2.0_RP,           0.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,         0.0_RP, -16.0_RP/11.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,         0.0_RP,          0.0_RP /),  &
  shape(RKCoef_a_7s6o_Butcher1964) )

 real(RP), parameter, public :: RKCoef_b_7s6o_Butcher1964(7) = &
     1.0_RP/120.0_RP * (/ 11.0_RP, 0.0_RP, 81.0_RP, 81.0_RP, -32.0_RP, -32.0_RP, 11.0_RP /)


  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------

  subroutine ATMOS_DYN_Tinteg_RKUtil_setup( this, rk_register_num, varname_list, &
      comm_id_offset, alloc_rkwork_flag )
    use scale_const, only: &
      UNDEF  => CONST_UNDEF
    use scale_comm_cartesC, only: &
      COMM_vars8_init      
    implicit none
    
    type(RKUtil), intent(inout) :: this
    integer, intent(in) :: rk_register_num
    character(*), intent(in) :: varname_list(:)
    integer, intent(in) :: comm_id_offset
    logical, intent(in) :: alloc_rkwork_flag

    integer :: reg_id
    integer :: var_id
    character(H_MID) :: rktag
    !------------------------------------------

    this%var_num = size(varname_list)
    this%register_num = rk_register_num

    if ( alloc_rkwork_flag ) call ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc( this )

    !-
    allocate( this%buf(KA,IA,JA,this%var_num) )
    allocate( this%comm_ind(this%var_num) )
    
    !$omp parallel workshare
    this%buf(:,:,:,:)   = UNDEF
    !$omp end parallel workshare

    do var_id = 1, this%var_num
      write(rktag,'(a,a,i2.2)') trim(varname_list(var_id)), 'RK'
      this%comm_ind(var_id) = comm_id_offset + var_id
      call COMM_vars8_init( trim(rktag), this%buf(:,:,:,var_id), this%comm_ind(var_id) )
    end do

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_setup

  subroutine ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc( this )
    implicit none

    type(RKUtil), intent(inout) :: this
    !------------------------------------------

    if (this%register_num > 0) then

      allocate( this%work0(KA,IA,JA,this%var_num) )
      allocate( this%work(KA,IA,JA,this%var_num,this%register_num) )

#ifdef DEBUG      
      !$omp parallel workshare
      this%work0 (:,:,:,:)   = UNDEF
      this%work  (:,:,:,:,:) = UNDEF
      !$omp end parallel workshare
#endif
    end if    

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc

  subroutine ATMOS_DYN_Tinteg_RKUtil_rkwork_dealloc( this )
    implicit none

    type(RKUtil), intent(inout) :: this
    !------------------------------------------

    if ( allocated(this%work) ) then
      deallocate( this%work0 )
      deallocate( this%work )
    end if    

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_rkwork_dealloc

  subroutine ATMOS_DYN_Tinteg_RKUtil_comm( this )
    use scale_comm_cartesC, only: COMM_vars8  
    implicit none

    type(RKUtil), intent(inout) :: this

    integer :: var_id
    !------------------------------------------

    do var_id = 1, this%var_num
      call COMM_vars8( this%buf(:,:,:,var_id), this%comm_ind(var_id) )
    end do 

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_comm

  subroutine ATMOS_DYN_Tinteg_RKUtil_comm_wait( this )
    use scale_comm_cartesC, only: COMM_wait
    implicit none

    type(RKUtil), intent(inout) :: this

    integer :: var_id
    !------------------------------------------

    do var_id = 1, this%var_num
      call COMM_wait( this%buf(:,:,:,var_id), this%comm_ind(var_id), .false. )
    end do 
    
    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_comm_wait

  !-------------------------
end module scale_atmos_dyn_tinteg_rkutil