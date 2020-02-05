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
    integer, allocatable :: comm_ind(:,:)
    real(RP), allocatable :: tmp(:,:,:,:,:)
    integer :: var_num
  end type

  public :: ATMOS_DYN_Tinteg_RKUtil_setup
  public :: ATMOS_DYN_Tinteg_RKUtil_comm
  public :: ATMOS_DYN_Tinteg_RKUtil_comm_wait

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

  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------

  subroutine ATMOS_DYN_Tinteg_RKUtil_setup( this, rk_register_num, varname_list, comm_id_offset )
    use scale_const, only: &
      UNDEF  => CONST_UNDEF
    implicit none
    
    type(RKUtil), intent(inout) :: this
    integer, intent(in) :: rk_register_num
    character(H_MID), intent(in) :: varname_list(:)
    integer, intent(in) :: comm_id_offset

    integer :: reg_id
    integer :: var_id
    character(H_MID) :: rktag
    !------------------------------------------

    this%var_num = size(varname_list)

    if (rk_register_num > 0) then
      allocate( this%tmp(KA,IA,JA,this%var_num,rk_register_num) )
      allocate( this%comm_ind(this%var_num,rk_register_num) )

      do reg_id = 1, rk_register_num
        do var_id = 1, this%var_num
          write(rktag,'(a,a,i2.2)') trim(varname_list(var_id)), 'RK', reg_id
          this%comm_ind(reg_id,var_id) = comm_id_offset + var_id
          call COMM_vars8_init( trim(rktag), this%tmp(:,:,:,var_id,reg_id), this%comm_ind(var_id,reg_id) )
        end do
      end do

      do reg_id = 1, rk_register_num 
        do var_id = 1, this%var_num
          this%tmp(:,:,:,var_id,reg_id) = UNDEF
        end do
      end do
    end if

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_setup

  subroutine ATMOS_DYN_Tinteg_RKUtil_comm( this, reg_id )
    implicit none
    type(RKUtil), intent(inout) :: this
    integer, intent(in) :: reg_id

    integer :: var_id
    !------------------------------------------

    do var_id = 1, this%var_num
      call COMM_vars8( this%tmp(:,:,:,reg_id,var_id), this%comm_ind(reg_id,var_id) )
    end do 

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_comm

  subroutine ATMOS_DYN_Tinteg_RKUtil_comm_wait( this, reg_id )
    implicit none
    type(RKUtil), intent(inout) :: this
    integer, intent(in) :: reg_id

    integer :: var_id
    !------------------------------------------

    do var_id = 1, this%var_num
      call COMM_wait( this%tmp(:,:,:,reg_id,var_id), this%comm_ind(reg_id,var_id), .false. )
    end do 
    
    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_comm_wait

  !-------------------------


end module scale_atmos_dyn_tinteg_rkutil